
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

J.test(){
echo "chr1	100	200	n	0	+	100	200	0	2	10,20	0,80
chr1	101	201	n	0	+	101	201	0	2	9,21	0,79
chr1	100	200	n	0	+	100	200	0	1	100	0
" | bed12_to_eej - 	
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
