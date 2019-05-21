
sam_to_bed12(){
usage="
$FUNCNAME [options] <bam> [region]
"; if [ $# -lt 1 ];then echo "$usage";return; fi
        samtools view -hb "$@" | bedtools bamtobed -bed12 -i stdin
}


count_junction(){
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
count_junction__test(){
echo "chr1	100	1000	g1	0	+	100	1000	0	3	100,200,300	0,200,600" \
| count_junction -
}

count_fragment(){
usage="$FUNCNAME <bed>"
if [ $# -lt 1 ];then echo "$usage $@"; return; fi
cat $1 | perl -e 'use strict;
	my %r=();
	while(<STDIN>){chomp;my@d=split/\t/,$_;
		$r{$d[0]}{$d[1]} |= 1;	
		$r{$d[0]}{$d[2]} |= 2;	
	}
	foreach my $c (keys %r){
		my @x=sort {$a<=>$b} keys %{$r{$c}};
		map{
			print $x[$_-1]," ",$x[$_],"\n";
		} 1..$#x;	
	}
'

}
count_fragment__test(){
echo "chr1	100	200	g1	1	+
chr1	300	400	g1	1	+
chr1	150	350	g1:intron	1	+" \
| count_fragment -	
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
