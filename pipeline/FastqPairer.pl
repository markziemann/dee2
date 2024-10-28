#!/usr/bin/env perl 

my $usage  = "Usage: FastqPairer.pl [-min <minLen>] [-symbol <symbol>] <FASTQ1> <FASTQ2>\n";

my $symbol;
my $minLen=0;

#Retrieve parameter
my @arg_idx=(0..@ARGV-1);
for my $i (0..@ARGV-1) {
	if ($ARGV[$i] eq '-symbol') {
		$symbol=$ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}elsif ($ARGV[$i] eq '-min') {
		$minLen=$ARGV[$i+1];
		delete @arg_idx[$i,$i+1];
	}
}
my @new_arg;
for (@arg_idx) { push(@new_arg,$ARGV[$_]) if defined $_; }
@ARGV=@new_arg;

my $fq1filename = shift or die $usage;
my $fq2filename = shift or die $usage;

# read the two files for getting paired read IDs
my %idCnt = ();

open(FILE,"<$fq1filename");
my $cnt=0;
while(<FILE>){
    $cnt++;
    
    if(($cnt%4)==1){
        chomp;
        $name=$_;
        if((defined $symbol) && /(.+)$symbol/){
            $name = $1;
        }elsif(/(.+)\s/){
            $name = $1;
        }
        $idCnt{$name}++;
    }
    if(($cnt%4)==2){
        chomp;
        if(length($_)<$minLen){
            $idCnt{$name}--;
        }
    }
}
close FILE;

open(FILE,"<$fq2filename");
$cnt=0;
while(<FILE>){
    $cnt++;
    
    if(($cnt%4)==1){
        chomp;
        $name=$_;
        if((defined $symbol) && /(.+)$symbol/){
            $name = $1;
        }elsif(/(.+)\s/){
            $name = $1;
        }
        $idCnt{$name}++;
    }
    if(($cnt%4)==2){
        chomp;
        if(length($_)<$minLen){
            $idCnt{$name}--;
        }
    }
}
close FILE;

# read and write for paired reads
open(FILE,"<$fq1filename");
open(OUT ,">$fq1filename.paired");
$cnt=0;
while(<FILE>){
    $cnt++;
    
    if(($cnt%4)==1){
        $line = $_;
        chomp $line;
        my $name=$line;
        if((defined $symbol) && $line=~/(.+)$symbol/){
            $name = $1;
        }elsif($line=~/(.+)\s/){
            $name = $1;
        }
        
        $flag=0;
        $flag=1 if $idCnt{$name}==2;
    }
    print OUT if $flag;
}
close FILE;
close OUT;

open(FILE,"<$fq2filename");
open(OUT ,">$fq2filename.paired");
$cnt=0;
while(<FILE>){
    $cnt++;
    
    if(($cnt%4)==1){
        $line = $_;
        chomp $line;
        my $name=$line;
        if((defined $symbol) && $line=~/(.+)$symbol/){
            $name = $1;
        }elsif($line=~/(.+)\s/){
            $name = $1;
        }
        
        $flag=0;
        $flag=1 if $idCnt{$name}==2;
    }
    print OUT if $flag;
}
close FILE;
close OUT;
