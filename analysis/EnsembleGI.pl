#! usr/bin/perl -w
use strict;

##This file will take a directory of Ensemble gene id files are parse them to make a single file.
#This file is a table containing the Ensemble Gene ID, Chromosome, Gene Start, Gene End, Complement (yes or no)

my $dir = "/Users/Chris/Google\ Drive/Uchicago/Nicolae/UsefulDocuments/EnsembleGeneAnnotations";
opendir(INPUT, $dir) || die "Cannot read directory $dir: $!";
my @dir_files = grep { (!/^\./) && -f "$dir/$_"} readdir INPUT;

my $output = "/Users/Chris/Desktop/Uchicago/Nicolae/GitWork/Methylation_DNAVariation/analysis/EnsembleGI.txt";

#Initialize Hash that stores all information#
#The keys are gene id and output is the vector of information in text format

my %Table;

##Read files and parse information##

foreach my $file (@dir_files) {
	open (my $fh, "< $dir/$file");
	my $chrom;							#Chromosom id
	my $gene;							#Gene id (also the key to hash)
	my $start;							#Gene start
	my $end;							#Gene end
	my $gene_cont = 0;					#Have I seen the first gene?
	my $comp = 0;						#Is it a complement gene?
	my $coding = 0;						#Is the gene a coding gene?
	my $gene_hit = 0;					#1 if next line is the gene id
	my $line;
	while ($line = <$fh>) {
		chomp $line;
		
		if ($line =~ /^ID/) {			#Initialize chromosome
			$chrom = Extract_Chrom(($line));
		}
		
		if (!$coding && $gene_cont && !$gene_hit && $line =~ /^FT[ ]{3}mRNA/) {
			$coding = 1;
		}
		
		if ($gene_hit) {				#Get gene name (also key id)
			$gene = Extract_Gene(($line));
			$gene_hit = 0;
		}		
		
		if ($line =~ /^FT[ ]{3}gene/) {   #Gene start and end
			if ($gene_cont) {
				if (!defined($Table{$gene}[0])) {
					@{$Table{$gene}} = ($chrom, $start, $end, $comp, $coding);
				} else {
					@{$Table{"$gene.$start.$end"}} = ($chrom, $start, $end, $comp, $coding);
				}		
			}
			
			($comp, $start, $end) = Extract_StartEnd(($line));
			$gene_cont = 1;
			$coding = 0;
			$gene_hit = 1;
		}
	}

	if (!defined($Table{$gene}[0])) {     #Record information for last gene
		@{$Table{$gene}} = ($chrom, $start, $end, $comp, $coding);
	} else {
		@{$Table{"$gene.$start.$end"}} = ($chrom, $start, $end, $comp, $coding);
	}
	
}


##Write %Table to a table##

open (OUT, "> $output") || die;
print OUT "EnsembleGI\tChromosome\tGene_Start\tGene_End\tComplement\tCoding\n";
foreach my $key (sort(keys(%Table))) {
	my $line_out = join("\t", @{$Table{$key}});
	print OUT "$key\t$line_out\n";
}
close OUT;



#######		Subroutines		#######


sub Extract_Chrom {
	my $line_sub = shift(@_);
	$line_sub =~ s/^ID[ ]+([A-Z,0-9]+)[ ]+.*$/$1/;
	return($line_sub);
}

sub Extract_StartEnd {
	my $line_sub = shift(@_);
	$line_sub =~ s/^FT[ ]+gene[ ]+//;
	my $comp_sub = 0;
	if ($line_sub =~ /complement/) {
		$comp_sub = 1;
		$line_sub =~ s/^complement\((.*)\)$/$1/;
	}
	my @out = split(/\.\./, $line_sub);
	return(($comp_sub, @out));
}

sub Extract_Gene {
	my $line_sub = shift(@_);
	$line_sub =~ s/^FT[ ]+\/gene=(.*)$/$1/;
	$line_sub =~ s/"//;
	return($line_sub);
}

#sub Extract_mRNA {
#	my $line_sub = shift(@_);
#	$line_sub =~ s/^FT[ ]+mRNA[ ]+(.*)$/$1/;
#	$line_sub =~ s/complement//;
#	$line_sub =~ s/join//;
#	$line_sub =~ s/\(//;
#	$line_sub =~ s/\)//;
#	return($line_sub);
#}
#
#sub Extract_mRNA_cont {
#	my $line_sub = shift(@_);
#	$line_sub =~ s/^FT[ ]+(.*)$/$1/;
#	$line_sub =~ s/complement//;
#	$line_sub =~ s/join//;
#	$line_sub =~ s/\(//;
#	$line_sub =~ s/\)//;
#	return($line_sub);
#}