
util.groupby(){
usage="
$FUNCNAME [options] <file> <grps> <cols> <opts>
 <grp>: columns used for grouping 
 <col> : column to be summarized 
 <opt> : operations
"

local OPTIND; local OPTARG;
local HEADER=0;
while getopts ":h" arg;do
case "$arg" in
	h) HEADER=1;;
	*) echo "$usage";return;;
esac
done 
shift $(( OPTIND - 1));

if [ $# -lt 4 ];then echo "$usage"; return; fi

cat $1 | perl -e 'use strict; 
	my $HEADER='$HEADER';
	my @G=(); my @C=(); my @OPT=();
	my $D=","; ## separator
	my %res=();
	my $first=1;
	my $NF; ## number of fields
	if( $HEADER ){my $line=<STDIN>;print "$line";}
	while(<STDIN>){chomp;my @a=split/\t/,$_;
		if($first){$first=0; $NF=scalar @a; 
			foreach my $e (split/,/,$ARGV[1]){
				if($e =~/(\d+)-(\d+)/){ push @G, $1..$2; }
				elsif($e =~/(\d+)-/){ push @G, $1..$NF;}
				elsif($e =~/-(\d+)/){ push @G, 0..$1;}
				else{ push @G,$e; }
			}
			my @col0=split/,/,$ARGV[2];
			my @opt0=split/,/,$ARGV[3];
			foreach my $i (0..$#col0){
				if($col0[$i] =~/(\d+)-(\d+)/){ push @C, $1..$2; map{ push @OPT,$opt0[$i] } $1..$2;}
				elsif($col0[$i] =~/(\d+)-/){ 
					push @C, $1..$NF; map{ push @OPT,$opt0[$i] } $1..$NF;}
				elsif($col0[$i] =~/-(\d+)/){ push @C, 0..$1; map{ push @OPT,$opt0[$i] } 0..$1;  }
				else{ push @C, $col0[$i]; push @OPT,$opt0[$i];} 
			}
			#print join(",",@G),"\n"; print join(",",@C),"\n"; print join(",",@OPT),"\n";
		}
		my $k=join("\t", map {$a[ $_ - 1 ]} @G);
		my @v= map {$a[ $_ - 1 ]} @C;
		if( defined $res{$k}){
			foreach my $i ( 0..$#C){
				if($OPT[$i] eq "collapse"){
					$res{$k}->[$i] .= $D.$v[$i];
				}elsif($OPT[$i] =~ /uniq/){
					my %kk=(); map{ $kk{$_}=1 } split /$D/,$res{$k}->[$i]; $kk{$v[$i]}=1;
					$res{$k}->[$i] = join($D,keys %kk);
				}elsif($OPT[$i] =~ /count/){
					$res{$k}->[$i] ++;
				}elsif($OPT[$i] eq "sum"){
					$res{$k}->[$i] += $v[$i];
				}elsif($OPT[$i] eq "min"){
					if( ! defined $res{$k}->[$i] || $res{$k}->[$i] > $v[$i] ){
						$res{$k}->[$i] = $v[$i];
					}
				}elsif($OPT[$i] eq "max"){ 
					if( ! defined $res{$k}->[$i] || $res{$k}->[$i] < $v[$i] ){
						$res{$k}->[$i] = $v[$i];
					}
				}else{
					die "unknown operation, \"$OPT[$i]\"";
				}
			}
		}else{
			foreach my $i ( 0..$#C){
				if($OPT[$i] eq "count"){ $v[$i] = 1; }
			}
			$res{$k}=\@v;
		}
	}
	foreach my $k (keys %res){
		print $k,"\t",join("\t",@{$res{$k}}),"\n";
	}
	
'  $@

}
util.groupby__test(){
echo \
"A	B	C	D	E
a	b	1	2	1
b	b	3	4	2
a	b	2	2	3
b	b	4	5	4" | util.groupby -h - 1,2 3- count

}
