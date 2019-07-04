simstudy(){
echo "
chr1	1000	2000	gene1	0	+	1000	2000	0,0,0	4	100,100,100,100	0,200,400,900
chr1	1000	2000	gene1	0	+	1000	2000	0,0,0	3	100,100,100	0,400,900
chr1	1000	1300	gene1	0	+	1000	1300	0,0,0	2	100,100	0,200
" | grep -v "^$" > tmp.transcript
salsa sim tmp.transcript 50 100 > tmp.read

}
