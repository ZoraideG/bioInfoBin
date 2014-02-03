use Bio::SearchIO; 
  use Bio::SeqIO;
 +use Getopt::Long;
  
  #basic clean up
  #sort and split genome
  
  
 -my $usage= "USAGE: extract.align.pl <genome> <blast.output> <TE_query.fas> <buffer> <seq number> <align/unalign>\n";
 -my $genome_in = shift or die $usage;
 -my $blast_in = shift or die $usage;
 -my $te_in = shift or die $usage;
 -my $buffer_in = shift or die $usage;
 -my $seq_num_in = shift or die $usage;
 -my $align = shift or die $usage;
 +
 +
 +my $genome_in 	= '';
 +my $blast_in 	= '';
 +my $te_in 		= '';
 +my $buffer_in 	= '';
 +my $seq_num_in 	= '';
 +my $align 		= ''; #sets align default to false
 +my $help		= '';
 +
 +GetOptions('genome=s' => \$genome_in,
 +		'blast=s' =>\ $blast_in,
 +		'consTEs=s' => \$te_in,
 +		'seqBuffer=i' => \$buffer_in,
 +		'seqNum=i' => \$seq_num_in,
 +        'align' => \$align, 
 +		'h'		=> \$help );
 +
 +if(@ARGV < 1){
 +	usage();
 +}
 +
 +	
 +
  
  
  #step 1 organizing blast output
 @@ -152,7 +170,28 @@ sub organize_blast_hits
  
  
  
 -
 +##################################################################################
 +sub usage
 +{
 +    print "\n";
 +	print "#####################################################################\n";
 +	print "#                                                                   #\n";	
 +	print "#  {{extract.align.pl - a program to extract TEs from a genome}}    #\n";
 +	print "#                                                                   #\n";
 +	print "#  Required values:                                                 #\n";                                              
 +	print "#    --genome      input genome file (in fasta format)              #\n";
 +	print "#    --blastInput  input blast file (tab format)                    #\n"; 
 +	print "#    --consTEs     file of consensus elements                       #\n";
 +	print "#    --seqBuffer   additional 5' and 3' bases in extracted sequence #\n";
 +	print "#    --seqNum      number of sequences to extract                   #\n";
 +	print "#                                                                   #\n";	
 +	print "#  Optional values:                                                 #\n";
 +	print "#    --align       aligns the sequences with MUSCLE                 #\n";
 +	print "#                                                                   #\n";
 +	print "#####################################################################\n";
 +	
 +	exit;
 +}	
  
  #my $size = (-s "$genome") / (1024 * 1024);
  #print $size;
