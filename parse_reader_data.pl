#Program to parse stacker-reader output files in .txt format and generate tab separated files by plate.
#!/usr/bin/perl
use strict;
use warnings;

opendir(DIR, "/BugScreen/") || die $!;
my @dir = readdir(DIR);

foreach my $f1(@dir){

	next unless $f1 =~ /(^NT.*)/;	#to read all bug folders starting with NT

	my $name1 = $1;

	opendir(NT, "/BugScreen/$f1") || die $!;
	my @contents = readdir(NT);

	foreach my $f2(@contents){

		next unless $f2 =~ /(^Replicate\d).*/;	#to read all Replicate folders starting with Replicate

		my $name2 = $1;

		opendir(REP, "/BugScreen/$f1/$f2") || die $!;
		my @rep = readdir(REP);

		foreach my $f3(@rep){

		next unless $f3 =~ /(^plate.*).txt/;

		my $name3 = $1;

		my $file = "/BugScreen/$f1/$f2/$f3";
		open my $fh, $file or die $!;
	
		my $outfile = "$name1"."_"."$name2"."_"."$name3".".tab";
		open(OUT, ">$outfile");

		while(my $line = <$fh>){

			#chomp $line;

			if ($line =~ /^Time\tT/ .. $line =~ /^Results/) {
        		
        		$line =~ s/,/\./g;

        		if($line =~ /^Results/){

        			next;
        		
        		}elsif($line =~ /^0:00:00/){

        			next;
        		
        		}elsif($line =~ /^\s+/){

        			next;

        		}else{

        		print OUT "$line";

        		}
   			 } 
		}	
	}
}
	
}	

#Below code removes the Temperature column from all files and formats time column

opendir(DIR, "/BugScreen/Data_analysis") || die $!;
my @dd = readdir(DIR);

foreach my $file(@dd){

	next unless $file =~ /(^NT.*\.tab)/;

	my $name = $1;

	my $ff  = `cut -f 1,3- $file`; #default delimiter for cut is tab

	if($name =~ /MIC/){

		$ff = `cut -f 1,3-74 $file`; #leave out G and H columns, they are empty wells
	}

	$ff =~ s/:\d\d:\d\d//g;

	open(OUT, ">$name") || die $!;

	print OUT $ff;

}
