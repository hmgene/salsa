SALSA_BINARY=$SALSA_HOME/bin/$( uname -sm | tr " " "." );
bedtools(){ 
	$SALSA_BINARY/$FUNCNAME $@
}
samtools(){ 
	$SALSA_BINARY/$FUNCNAME $@
}
