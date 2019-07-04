

countunsplice(){
usage="
$FUNCNAME <F> <read.bed6>
"
if [ $# -lt 2 ];then echo "$usage";return; fi
        hm bed.intersect $1 $2  -wa -wb  \
        | perl -e 'use strict; my %r=(); my $ncol=-1;
                sub getv{ my ($h)=@_; return defined $h ? $h : 0;}
                while(<STDIN>){ chomp; my @d=split/\t/,$_;
                        my $k=join("\t",@d[0..5]);
                        if( $d[1] > $d[7] ){ $r{ $k }{ 0 } ++; }
                        if( $d[2] < $d[8] ){ $r{ $k }{ 2 } ++; }
                        $r{ $k }{ 1 } ++;
                }
                foreach my $k (sort keys %r){
                        print $k,"\t";
                        print join("\t", map{ getv( $r{ $k }{ $_ }) } 0..2),"\n";
                }
        '
}
countunsplice__test(){
echo "chr1	1000	2000	ex1	0	+" > tmp.a
echo "chr1	1000	1002	1
chr1	1000	1001	2
chr1	1999	2001	10
chr1	999	1001	3" > tmp.b
countunsplice tmp.a tmp.b
rm tmp.*

}

testevent(){
usage="$FUNCNAME <ctr>[,<ctr>] <trt>[,<trt>] <outd>"
if [ $# -lt 3 ];then echo "$usage";return;fi
	local ctr=`perl -e 'my@d=('$1'); print join(",",map{"ctr$_"} 0..$#d);'`;
	local trt=`perl -e 'my@d=('$2'); print join(",",map{"trt$_"} 0..$#d);'`;

        perl -e 'use strict;
        my @ctr=split/,/,"'$1'";
        my @trt=split/,/,"'$2'";
	my %r=();
	map { my $fid=$_;
		open(my $fh,"<",$ctr[$_]) or die "$!";	
		while(<$fh>){chomp; my @d=split/\t/,$_; $r{join("\t",@d[0..4])}{$fid}=$d[5]; }
		close($fh);
	} 0..$#ctr;

	map { my $fid=$_+$#ctr+1;
		open(my $fh,"<",$ctr[$_]) or die "$!";	
		while(<$fh>){chomp; my @d=split/\t/,$_; $r{join("\t",@d[0..4])}{$fid}=$d[5]; }
		close($fh);
	} 0..$#trt;
	my $eid=0;
	my $gid=0;
	my $fid=0;
	print "chrom\tgene\tstrand\ttype\tgene_id\tfeature_id\tfeature\t";
	print join("\t",map{ "ctr$_" } 0..$#ctr),"\t";
	print join("\t",map{ "trt$_" } 0..$#trt),"\n";
	foreach my $k (keys %r){
		my ($c,$g,$t,$s,$tmp)=split/\t/,$k;
		my @fea=split/,/,$tmp;
		foreach my $i (0..$#fea){	
			print "$c\t$g\t$s\t$t\tg$gid\tf$fid\t$fea[$i]";
			foreach my $j (0..($#ctr + $#trt + 1)){
				my @tmp=split/,/,$r{$k}{$j};
				my $v=defined $tmp[$i] ? $tmp[$i] : 0;
				print "\t$v";
				$fid++;
			}
			print "\n";
		}
		$gid++;
	}
        ' | hm drimseq.run - $ctr $trt $3
}
testevent__test(){
echo "chr1	g1	+	SE	100^200,100^400,300^400	1,2,3" > tmp.a
testevent tmp.a,tmp.a tmp.a tmp.o
rm -r tmp.*
}


nov(){
	local tmpd=`mktempd`;
	cat $1 > $tmpd/a
	salsa bedtools intersect -a $tmpd/a -b $tmpd/a -wao \
	| perl -e 'use strict; 
	my %r=();
	while(<STDIN>){chomp; my @d=split/\t/,$_;
		my $i=int($#d/2);
		my $k=join("\t",@d[0..($i-1)]);
		if(!defined $r{$k}){ $r{$k}=0;}
		if($d[$#d] < $d[2]-$d[1]){ $r{$k}++; }
	}
	print join("\n",grep {$r{$_} == 0 } keys %r),"\n";
	'
	rm -r $tmpd
}
nov__test(){
echo "chr1	1	200
chr1	100	200
chr1	59	200
chr1	50	60" | nov -
}
makeRIevent(){
usage="
$FUNCNAME <exon> <j> <bdg> 
"
if [ $# -lt 3 ];then echo "$usage";return; fi
{
	awk '{print "1\t"$0}' $1
	nov $2 | awk '{print "2\t"$0}'
	awk '{print "3\t"$0}' $3
} | perl -e 'use strict; 
	my %A=();
	my %B=();
	while(<STDIN>){chomp;my($f,@d)=split/\t/,$_;
		if($f==1){
			my $k=$d[0]."\t".$d[3]."\t".$d[5];
			$B{$k}{$d[1]}++;
			$B{$k}{$d[2]}++;
		}elsif($f==2){
			$A{$d[0]}{$d[1]}{0}{$d[2]}+=$d[3];
			$A{$d[0]}{$d[2]}{1}{$d[1]}+=$d[3];
		}else{
			if(defined $A{$d[0]}{$d[1]}{0}){
				$A{$d[0]}{$d[1]}{2}=$d[2];
			} 
			if(defined $A{$d[0]}{$d[2]}{1}){
				$A{$d[0]}{$d[2]}{3}=$d[1];
			}
		}
	}
	foreach my $k (keys %B){
		my ($c,$g,$t)=split/\t/,$k;
		my @x=sort {$a<=>$b} keys %{$B{$k}};
		foreach my $i (0..($#x-1)){ 
		foreach my $j (($i+1)..$#x){
			my $ij=$A{$c}{$x[$i]}{0}{$x[$j]};
			my $ii= $A{$c}{$x[$i]}{2};
			my $jj= $A{$c}{$x[$j]}{3};
			if( defined $ij && defined $ii && defined $jj){
				print "$c\t$g\t$t\tRI\t$x[$i]^$x[$j],$x[$i]-$ii,$jj-$x[$j]\n";
			}
		}}
	}
'
}

makeRIevent__test(){
echo "chr1	1000	2000	g1	t1	+	1000	2000	0	2	100,50	0,950
chr1	1000	2000	g1	t2	+	1000	2000	0	1	1000	0" \
> tmp.gene
salsa flattenexon tmp.gene > tmp.exon
sim tmp.gene > tmp.read
countjunction tmp.read > tmp.junction
bed12toexon tmp.read | bdg - > tmp.bdg
makeRIevent tmp.exon tmp.junction tmp.bdg
rm tmp.*
}

countevent(){
usage="
$FUNCNAME <event> <j> <bdj>
"
if [ $# -lt 2 ];then echo "$usage";return; fi
perl -e 'use strict;
	sub getv{ my ($h)=@_; return defined $h ? $h : 0;}
	my %r=();
	my %D=();
	my @files=qw( '$1' '$2' '$3');
	my $fh; foreach my $file (@files){
		if($file eq "-"){ $fh=*STDIN; }else{ open($fh,"<",$file) or die "$!";}
		while(<$fh>){chomp; my @d=split/\t/,$_;
			if($file eq $files[0]){
				$r{join("\t",@d[0..4])}=1;
			}elsif($file eq $files[1]){
				$D{$d[0]}{ $d[1]."^".$d[2] } += $d[3];
			}elsif($file eq $files[2]){
				$D{$d[0]}{ $d[1]."-".$d[2] } += $d[3];
			}
		}
		close($fh) unless $file eq "-";
	}

	foreach my $k (keys %r){
		my ($c,$g,$t,$type,$tmp)=split/\t/,$k;
		my @j=split/,/,$tmp;
		print $k,"\t",join(",", map { getv( $D{$c}{$_} ) } @j),"\n";
	}
	'

	
}
countevent__test(){
echo "chr1	100	200	g1	00	+
chr1	200	300	g1	00	+
chr1	400	500	g1	00	+
chr1	600	700	g1	00	+
chr1	300	400	g1	00	+" > tmp.a

echo "chr1	200	300	1
chr1	200	400	2
chr1	500	600	4
chr1	200	600	3" > tmp.b

echo "chr1	200	210	1
chr1	280	300	1" > tmp.c
makeevent tmp.a tmp.b tmp.c > tmp.d
countevent tmp.d tmp.b tmp.c
rm tmp.*
}
makeevent(){
usage="
$FUNCNAME <exon> <j> <bdg>
"
if [ $# -lt 3 ];then echo "$usage";return; fi
	perl -e 'use strict;
		my @files=qw( '$1' '$2' '$3' );
		my %g=();
		my %j=();
		my %B=();
		my $fh; foreach my $file (@files){
			if($file eq "-"){ $fh=*STDIN; }else{ open($fh,"<",$file) or die "$!";}
			while(<$fh>){chomp; my @d=split/\t/,$_;
				if($file eq $files[0]){
					my $k=join("\t",$d[0],$d[3],$d[5]);
					$g{$k}{$d[1]}=1;
					$g{$k}{$d[2]}=1;
				}elsif($file eq $files[1]){
					$j{$d[0]}{$d[1]}{$d[2]} += $d[3];
					$j{$d[0]}{$d[2]}{$d[1]} += $d[3];
				}elsif($file eq $files[2]){
					$B{$d[0]}{$d[1]}{$d[2]} += $d[3];
					$B{$d[0]}{$d[2]}{$d[1]} += $d[3];
				}		
			}
			close($fh) unless $file eq "-";
		
		}
		foreach my $k (keys %g){
			my ($c,$n,$t)=split /\t/,$k;
			my @x=sort {$a<=>$b} keys %{$g{$k}};
			map{ my $w=$_;
				my @e=sort {$a<=>$b} grep {defined $g{$k}{$_} && $_ > $w } keys %{$j{$c}{$w}};
				if($#e > 0){
					my $type = $t eq "+" ? "A3" : "A5";
					print "$c\t$n\t$t\t$type\t",join(",",map{"$w^$_"} @e),"\n";	

					##   w]----[x    y]----[z   
					my $x=$e[0];	
					foreach my $z (@e[1..$#e]){
					foreach my $y ( grep { $_ > $x && $_ < $z} keys %{$j{$c}{$z}}){
						print "$k\tSE\t$w^$x,$w^$z,$y^$z\n";
					}}
				}

				
				@e=sort {$a<=>$b} grep {defined $g{$k}{$_} && $_ < $w } keys %{$j{$c}{$w}};
				if($#e > 0){
					my $type = $t eq "+" ? "A5" : "A3";
					print "$c\t$n\t$t\t$type\t",join(",",map{"$_^$w"} @e),"\n";	
				}
			} @x[1..($#x-1)];
			## RI
			foreach my $i (1..($#x-2)){ 
			foreach my $j (($i+1)..($#x-1)){ 
				my $hit=0;
				map{ my $ii=$_;
					map { my $jj=$_;
						print "$c\t$n\t$t\tRI\t$x[$i]^$x[$j],$x[$i]-$ii,$jj-$x[$j]\n";
						$hit=1;
					} grep {$_ < $x[$j] } keys %{$B{$c}{$x[$j]}};
				} grep {$_ > $x[$i]} keys %{$B{$c}{$x[$i]}};
				next if($hit); ## ignore nesting RI events
			}}
		}
	
	'
}

makeevent__test(){
echo "chr1	100	200	g1	00	+
chr1	200	300	g1	00	+
chr1	400	500	g1	00	+
chr1	600	700	g1	00	+
chr1	300	400	g1	00	+" > tmp.a

echo "chr1	200	300	1
chr1	200	400	2
chr1	500	600	4
chr1	200	600	3" > tmp.b
echo "chr1	200	210	1
chr1	280	300	1" > tmp.c
makeevent tmp.a tmp.b tmp.c
rm tmp.*
}


countexonfragment(){
usage="
$FUNCNAME <exonfragment> <bdg>
"
if [ $# -lt 2 ];then echo "$usage";return; fi
	bedtools intersect -a ${1/^-$/stdin} -b ${2/^-$/stdin} -wa -wb  \
	| perl -e 'use strict; my %r=(); my $ncol=-1;
		sub getv{ my ($h)=@_; return defined $h ? $h : 0;}
		while(<STDIN>){ chomp; my @d=split/\t/,$_;
			print $_,"\n";
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

flattenexon(){
usage="
$FUNCNAME <gene.bed12> <junction.bdg> [<junction.bdg>..]
"
if [ $# -lt 1 ];then echo "$usage";return; fi
	{
		bed12toexon $1
		if [ $# -gt 1 ];then
			cat ${@:2}
		fi
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
				## check directionality 
				if( $y > 0 && $_ < $x || $ y< 0 && $_ > $x ){ 
					$h{$_}=1; ## found novel (exon-intron) junction
				}
			} grep { ! defined $r{$k}{$_} } keys %{$j{$c}{$x}};
			
		} keys %h;
		my @z=sort {$a<=>$b} keys %h;
		map  {
			print join("\t",$c,$z[$_-1],$z[$_],$g,"$h{$z[$_-1]}$h{$z[$_]}",$t),"\n";
		} grep { #filter out beyond gene boundary
			!( $_==1 && $h{$z[$_-1]} == 1 || $_==$#z && $h{$z[$_]} == 1 )
		} 1..$#z;
	}
	'
}
flattenexon__test(){
## [100 200)---[300 500)
echo "chr1	100	1000	g1	0	+	100	1000	0	3	100,200,300	0,200,600" > tmp.a
echo "chr1	200	250	1
chr1	280	300	1
chr1	50	100	1
chr1	50	300	2
chr1	500	600	1
" > tmp.b
flattenexon tmp.a tmp.b
rm tmp.*
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
                } 0..($d[9]-1);
		
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

bdg(){
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
bdg__test(){
echo "chr1	100	200	g1	1	+
chr1	300	400	g1	1	+
chr1	150	250	g1:intron	1	+" \
| bdg -	
}

sim(){
usage="
$FUNCNAME [options] <transcript.bed12> 
 [optinos]:
	-l <readlen>
	-n <readnum>
"
local OPTIND; local OPTARG;
n=100; l=100
while getopts ":l:n:h" arg;do
echo "$l $n"
case "$arg" in
        l) l="${OPTARG}";;
        n) n="${OPTARG}";;
        *) echo "$usage";return;;
esac     
done  
shift $(( OPTIND - 1));

if [ $# -lt 1 ];then echo "$usage $@"; return; fi
	cat $1 | perl -e 'use strict;
		my $l_read='$l';
		my $n_read='$n';
		my %T=();
		my $total=0;
		while(<STDIN>){ chomp; my @d=split/\t/,$_;
			my $l=0;map { $l+=$_ } split /,/,$d[10];
			next if($_ eq "");
			$T{$_} = $l;
			$total += $l;
		}
		foreach my $t (keys %T){
			my $p=$T{$t}/$total*$n_read;
			my@d=split/\t/,$t;
			## get a transcript
			my @l=split/,/,$d[10]; 
			my $n=0;map{ $n+=$_;} @l;
			my @s=map { $d[1] + $_ } split/,/,$d[11];		

			## generate reads
			my @r_s=sort {$a<=>$b} map { int(rand($n-$l_read)) } 0..int($p);
			my @r_e=map { $_+$l_read } @r_s;

			## map to genomic positions
			my @t=@l; map { $t[$_]+=$t[$_-1]; } 1..$#t;
			my $i=0;
			my %h=();
			foreach my $x (sort { $a<=>$b} @r_s, map{ $_ -1 } @r_e ){
				while( $t[$i]<=$x && $i<=$#t){ $i++; }
				$h{ $x }{i} = $i;
				$h{ $x }{y} = $x - ($i> 0 ? $t[ $i-1] :0 ) + $s[ $i];
			}
			map {
				my ($i,$j)=map{ $h{ $_ }{i} } ($r_s[$_], $r_e[$_]-1);
				my ($start,$end)=map{ $h{ $_ }{y} } ($r_s[$_], $r_e[$_]-1);
				my @ss=();
				my @ee=();
				if($i==$j){
					push @ss,$start; 
					push @ee,$end+1;
				}else{
					push @ss,$start;
					push @ee,$s[$i]+$l[$i];
					map {
						push @ss,$s[$i];
						push @ee,$s[$i]+$l[$i];
					} ($i+1)..($j-1);
					push @ss,$s[$j];
					push @ee,$end+1;
				}
				my $sizes=join(",",map{ $ee[$_] - $ss[$_] } 0..$#ss);
				my $starts=join(",",map{ $ss[$_] - $ss[0] } 0..$#ss);
				print join("\t",$d[0],$ss[0],$ee[$#ee],@d[3..5],$ss[0],$ee[$#ee],0,$#ss+1,$sizes,$starts),"\n";
			} 0..$#r_s;
		}
		
	'
}
sim__test(){
echo "chr1	100	1000	g1	0	+	100	1000	0	3	100,200,300	0,200,600
chr1	100	1000	g2	0	+	100	1000	0	2	100,300	0,600" \
| salsa sim -
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
