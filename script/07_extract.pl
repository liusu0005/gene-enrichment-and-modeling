
use strict;
#use warnings;

my %hash=();

my %aliasHash=();

open(RF,"../output/06/intersect.txt") or die $!;
while(my $line=<RF>){
	chomp($line);
	if(exists $aliasHash{$line}){
		$hash{$aliasHash{$line}}=1;
	}
	else{
	  $hash{$line}=1;
  }
}
close(RF);

open(RF,"../sample_data/geneMatrix.txt") or die $!;
open(WF,">../output/07/output.txt") or die $!;
while(my $line=<RF>){
	if($.==1){
		print WF $line;
		next;
	}
	my @arr=split(/\t/,$line);
	my @zeroArr=split(/\|/,$arr[0]);
	if(exists $hash{$zeroArr[0]}){
		print WF $line;
		delete($hash{$zeroArr[0]});
	}
}
close(WF);
close(RF);

foreach my $key(keys %hash){
	print $key . "\n";
}

