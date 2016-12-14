salsa.junction(){
usage="
$FUNCNAME <bed12>
";if [ $# -lt 1 ];then echo "$usage"; return; fi
cat $1 | perl -e 'use strict; my %res=();
	while(<STDIN>){ chomp; my @a=split/\t/,$_;
		print $_,"\n"; 
		my @l=split/,/,$a[10];
		my @s=split/,/,$a[11];
		for(my $i=0; $i < $a[9] - 1; $i++){
			my $start= $a[1] + $s[$i] + $l[$i] - 1;
			my $end = $a[1] + $s[$i+1] + 1;
			$res{$a[5]}{$a[0]}{$start."\t".$end} ++;
		}
	}
	foreach my $st (keys %res){
	foreach my $c (keys %{$res{$st}}){
	foreach my $se (keys %{$res{$st}{$c}}){
		print join("\t",( $c,$se,"j",$res{$st}{$c}{$se},$st)),"\n";
	}}}
'
} 
salsa.junction.test(){
echo "
01234567890123456789
 EEEE-EEE---EEEEEE
" | salsa gen_bed12 - | salsa junction -
}

salsa.gen_bed12(){
usage="$FUNCNAME <intput.txt>"
if [ $# -lt 1 ];then echo "$usage"; return; fi
cat $1 | perl -e 'use strict; 
	sub f{
		my ($s,$e,$g) = @_;
		my $n=scalar @$s;
		#print join(",",@$s)," ",join(",",@$e),"\n";
		my $s0=$s->[0];		
		my $e0=$e->[$n-1];		
		my $st="+"; if( lc $g eq $g){ $st="-";}
		print "chr1\t$s0\t$e0\t$g\t0\t$st\t$s0\t$e0\t0,0,0\t$n\t";
		print join(",",map { $e->[$_] - $s->[$_]} (0..($n-1))),"\t";
		print join(",",map { $_ - $s0 } @$s),"\n";

	}
	while(<STDIN>){chomp;  my $tmp=" ".$_." "; $tmp=~s/\d/ /g;
		my @S=split//,$tmp;
		my @sizes=(); my @starts=(); my @ends=();
		my $type="";
		my $start=0;
		for(my $i=1; $i< scalar @S; $i++){ 
			my ($a,$b)=($S[$i-1],$S[$i]);
			if($a eq " " && $b =~/\w/){
				push @starts,$i-1;
			}elsif($a =~/\w/ && $b eq " "){
				push @ends,$i-1;
				f(\@starts,\@ends,$a);
				@starts=();@ends=();
			}elsif($a eq "-" && $b =~/\w/){
				push @starts,$i-1;
			}elsif($a =~/\w/ && $b eq "-"){
				push @ends,$i-1;
			}elsif($a ne $b){
				push @ends,$i-1;
				f(\@starts,\@ends,$a);
				@starts=();@ends=();
				push @starts,$i-1;
			} 
		}
	}
'
}



