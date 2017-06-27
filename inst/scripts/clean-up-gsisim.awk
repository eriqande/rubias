
# print column headers that we will need
BEGIN {
	printf("collection  indiv ");
	go = 0
}

# print headers for locus names
NF == 1 { printf(" %s %s", $1, $1); next}

# when you are starting a Pop block, get the pop name, etc
$1 == "POP" || $1 == "Pop" || $1 == "pop" {
	if(go == 0) {
		printf("\n");  # at first Pop print a return
	}
	go = 1;
	pop = $2;
	next;
}

# print the genotype info in a row
go == 1 && NF > 0 {
	printf("%s  ", pop);
	print $0
}

