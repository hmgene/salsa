countexonfragment(){
usage="
$FUNCNAME <exonfragment> <bdg>
"
if [ $# -lt 2 ];then echo "$usage";return; fi
	bedtools intersect -a ${1/^-$/stdin} -b ${2/^-$/stdin} -wa -wb  \
	| perl -e 'use strict; my %r=(); my $ncol=-1;
		sub getv{ my ($h)=@_; return defined $h ? $h : 0;}
		while(<STDIN>){ chomp; my @d=split/\t/,$_;
			my $k=join("\t",@d[0..5]);
			if( $d[1] > $d[7] ){ $r{ $k }{ 0 } +=$d[9]; }
			if( $d[2] < $d[8] ){ $r{ $k }{ 2 } +=$d[9]; }
			$r{ $k }{ 1 } += $d[9];
		}
		foreach my $k (sort keys %r){
			print $k,"\t";
			print join("\t", map{ getv( $r{$k}{$_} ) } 0..2),"\n";
		}
	'
}
countexonfragment__test(){
echo "chr1	100	200	g1	00	+
chr1	200	300	g1	00	+
chr1	300	400	g1	00	+" > tmp.a

echo "chr1	199	200	1
chr1	199	201	2
chr1	300	400	3" > tmp.b
countexonfragment tmp.a tmp.b
rm tmp.*

}

makeexonfragment(){
usage="
$FUNCNAME <gene.bed12> <junction.bdg> [<junction.bdg>..]
"
if [ $# -lt 1 ];then echo "$usage";return; fi
	{
		bed12toexon $1
		cat ${@:2}
	} | perl -e 'use strict; 
	my %r=();
	my %j=();
	while(<STDIN>){chomp; my @d=split/\t/,$_;
		if($#d == 5){ ## exon
			$r{join("\t",$d[0],$d[3],$d[5])}{$d[1]}++; 
			$r{join("\t",$d[0],$d[3],$d[5])}{$d[2]}--; 
		}elsif($#d==3){ ## junction
			$j{$d[0]}{$d[1]}{$d[2]} += $d[3];
			$j{$d[0]}{$d[2]}{$d[1]} -= $d[3];
		}
	}
	## produce overlapping fragments of different genes
	foreach my $k (keys %r){
		my ($c,$g,$t)=split/\t/,$k;
		my %h=map {$_=>0} keys %{$r{$k}}; ## exon boundaries
		map { my ($x,$y)=($_,$r{$k}{$_});
			map { ## $_^[x  ] or  [  x]^$_
				if( $y > 0 && $_ < $x || $ y< 0 && $_ > $x ){
					$h{$_}=1; ## found novel (exon-intron) junction
				}
			} grep { ! defined $r{$k}{$_} } keys %{$j{$c}{$x}};
			
		} keys %h;
		my @z=sort {$a<=>$b} keys %h;
		map  {
			print join("\t",$c,$z[$_-1],$z[$_],$g,"$h{$z[$_-1]}$h{$z[$_]}",$t),"\n";
		} 1..$#z;
	}
	'
}
makeexonfragment__test(){
## [100 200)---[300 500)
echo "chr1	100	1000	g1	0	+	100	1000	0	3	100,200,300	0,200,600" > tmp.a
echo "chr1	200	250	1
chr1	280	300	1
" > tmp.b
makeexonfragment tmp.a tmp.b
rm tmp.*
}

makeevent(){
usage="
$FUNCNAME <gene.bed12> <junction.bdg> [..<junction.bdg] 
"; if [ $# -lt 2 ];then echo "$usage";return; fi
        samtools view -hb "$@" | bedtools bamtobed -bed12 -i stdin
}


samtobed12(){
usage="
$FUNCNAME [options] <bam> [region]
"; if [ $# -lt 1 ];then echo "$usage";return; fi
        samtools view -hb "$@" | bedtools bamtobed -bed12 -i stdin
}
bed12toexon(){
usage="
$FUNCNAME <bed12> 
"; if [ $# -lt 1 ];then echo "$usage";return; fi
	cat $1 | perl -ne 'chomp; my @d=split/\t/,$_;
                my @l=split/,/,$d[10]; my @s=split/,/,$d[11];
                map {
			print join("\t",$d[0],$d[1]+$s[$_],$d[1]+$s[$_]+$l[$_],@d[3..5]),"\n";
                } 0..($d[9]-2);
		
	'
}
bed12toexon__test(){
echo "chr1	100	1000	g1	0	+	100	1000	0	3	100,200,300	0,200,600" \
| bed12toexon -
}


countjunction(){
usage="
$FUNCNAME <bed12>
"; if [ $# -lt 1 ];then echo "$usage";return; fi
   cat $1 | perl -e 'my %r=();
        while(<STDIN>){chomp;my @a=split/\t/,$_;
                my @l=split/,/,$a[10]; my @s=split/,/,$a[11];
                map {
                        $r{ $a[0] }{ $a[1]+$s[$_]+$l[$_] }{ $a[1]+$s[$_+1] } ++;
                } 0..($a[9]-2);
        }
        foreach my $c (keys %r){
	foreach my $s (keys %{$r{$c}}){
	foreach my $e (keys %{$r{$c}{$s}}){
                print join("\t",$c,$s,$e,$r{$c}{$s}{$e}),"\n";
        }}}
        '

}
countjunction__test(){
echo "chr1	100	1000	g1	0	+	100	1000	0	3	100,200,300	0,200,600" \
| countjunction -
}

countfragment(){
usage="$FUNCNAME <bed>"
if [ $# -lt 1 ];then echo "$usage $@"; return; fi
cat $1 | perl -e 'use strict;
	my %r=();
	while(<STDIN>){chomp;my@d=split/\t/,$_;
		$r{$d[0]}{$d[1]}++;
		$r{$d[0]}{$d[2]}--;
	}
	foreach my $c (keys %r){
		my @x=sort {$a<=>$b} keys %{$r{$c}};
		my $acc=0;
		map{ my $e=$x[$_]; my $s=$x[$_-1];
			$acc+= $r{$c}{$s};			
			if($acc > 0){
				print "$c\t$s\t$e\t$acc\n";
			}
		} 1..$#x;	
	}
'

}
countfragment__test(){
echo "chr1	100	200	g1	1	+
chr1	300	400	g1	1	+
chr1	150	250	g1:intron	1	+" \
| countfragment -	
}

sim(){
usage="$FUNCNAME <bed12>"
if [ $# -lt 1 ];then echo "$usage $@"; return; fi
	cat $1 | perl -e 'use strict;
		my $l_read=50;
		my $n_read=100;
		while(<STDIN>){ chomp; my@d=split/\t/,$_;
			## get a transcript
			my @l=split/,/,$d[10]; 
			my @ll=(0,@l[0..($#l-1)]);
			map{$ll[$_]+=$ll[$_-1];} 1..$#ll;
			my @s=split/,/,$d[11];		

			## generate reads
			my @r_s=sort {$a<=>$b} map { int(rand($ll[$#ll]-$l_read)) } 0..$n_read;
			my @r_e=map { $_+$l_read } @r_s;

			## map to genomic positions
			my %h=(); my $i=0; map{
				$h{ $_ } = $i; ## exon number 
				print "$i $_\n";
				if($_ >= $ll[$i] ){ $i++;}
			} sort {$a<=>$b} (@s, @r_s, @r_e);

			map {	
				my ($s1,$e1)=($r_s[$_],$r_e[$_]);
				map{ 
					my $s2=$d[1]+$s[$_]+($s1-$ll[$_]);
					#print "$_:$ll[$_] $s1 -> $s2\n";
				} $h{$s1};
			} 0..$#r_s;

		}
		
	'
}
sim__test(){
echo "chr1	100	1000	g1	0	+	100	1000	0	3	100,200,300	0,200,600" \
| sim -
}

R(){
usage="
$FUNCNAME <sam|bam>
"
if [ $# -lt 1 ];then echo "$usage"; return; fi
        samtools view -hb $@ | bedtools bamtobed -bed12 -i stdin
}

J(){ usage="
$FUNCNAME <sam|bam> 
"; if [ $# -lt 1 ];then echo "$usage"; return; fi
	R $1 | bed12_to_eej > $2.eej
	bedtools intersect -a $2.eij -b $1 -wa -wb
}

unspliced(){
usage="
$FUNCNAME <splice_junction> <read.bed12>
"; if [ $# -ne 2 ];then echo "$usage"; return; fi

		
}

bed12_to_eej(){
usage="
$FUNCNAME <read.bed12>
"
if [ $# -ne 1 ];then echo "$usage"; return; fi
	perl -e 'my %r=();
	while(<STDIN>){chomp;my @a=split/\t/,$_;
		my @l=split/,/,$a[10];
		my @s=split/,/,$a[11];
		map { 
			my @b=map { $a[1] + $_ } ( $s[$_]+$l[$_], $s[$_+1] );
			$r{ join("\t",$a[0],$a[1]+$s[$_]+$l[$_],$a[1]+$s[$_+1]) }{ $l[$_]."&".$l[$_+1] } ++ ;
		} 0..($a[9]-2);
	}
	foreach my $k (keys %r){
		print $k,"\t";
		print join(",", map { $_."=".$r{$k}{$_} } keys %{$r{$k}} ),"\n";
	}
	'
}

bed12_to_eej__test(){
echo "chr1	100	200	n	0	+	100	200	0	2	10,20	0,80
chr1	101	201	n	0	+	101	201	0	2	9,21	0,79
chr1	100	200	n	0	+	100	200	0	1	100	0
" | bed12_to_eej -	
}
