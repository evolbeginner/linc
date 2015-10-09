#! /usr/bin/perl

package get_WGD_block;

require Exporter;
use strict;
use 5.010;

our @ISA = qw(Exporter);
our @EXPORT = qw(get_WGD_block_Bowers_2003 get_WGD_block_regular);


sub get_WGD_block_regular{
	my ($infile, $filter_S, $type, $numOfGenesInBlockCutoff) = @_;
	open(my $IN, '<' , "$infile") or die "WGD_block_file cannot be opened";
	my (%strand, %block_locus, %block_locus_order);
    my (%block_counter);
	
	while(my $line=<$IN>){
		chomp($line);
		my (%symbol, $block_locus_core, $block, $symbol1, $symbol2);
		my ($block_full_locus_name, $block);
		# At1g02000	At4g00110	A01N036	(可能有一个或多个其他field)
	
        if ($type eq 'PGDD'){
            ($block, $symbol1, $symbol2) = (split /,/, $line)[0,3,4];
        }
        else{
            ($block_full_locus_name, $symbol1, $symbol2) = (split /\s+/, $line)[0,1,2];
        }
       
        
        #($symbol1,$symbol2) = sort ($symbol1,$symbol2);
        ($symbol1,$symbol2) = ($symbol1,$symbol2);
 
		for ($symbol1,$symbol2){
			$_ = uc($_);
		}

		given ($type){
			when($type eq 'Bowers'){
				$block_locus_core = $block_full_locus_name;
				($block) = $block_locus_core =~ /^(.\d\d)/;
			}
            when($type eq 'PGDD'){
                my $counter = ++$block_counter{$block};
                $block_locus_core = $block . '-' . $counter;
            }
			default{
				my $order = $block_locus_order{$block_full_locus_name}++;
				$block_locus_core = $block_full_locus_name.'.'.$order;
			}
		}
        next if not $block;
		next if exists $block_locus{$block}{$block_locus_core};
		#next if $block_full_locus_name =~ '^A\d+OA';
		$block_locus{$block}{$block_locus_core} = $symbol1.'!'.$symbol2;
	}

    foreach my $block (keys %block_counter){
        if($block_counter{$block} < $numOfGenesInBlockCutoff){
            delete $block_locus{$block};
        }
    }

	return(\%block_locus);
}


################################################################################
#          *************************** 废物 **************************         #
sub get_WGD_block_PGDD{
	my ($infile, $filter_S, $type) = @_;
	open(my $IN, '<' , "$infile") or die "WGD_block_file cannot be opened";
	my (%strand, %block_locus, %block_locus_order);
	
	while(<$IN>){
        # BLOCK_NO,BLOCK_SCORE,E_VALUE,LOCUS_1,LOCUS_2,Ka,Ks
        # 231,236.0,0.0,16052055,ATCG00860,0.0026,0.007
		chomp;
		my (%symbol, $block_locus_core, $block, $symbol1, $symbol2);
		my ($block_full_locus_name);
		($block_full_locus_name, $symbol1, $symbol2) = (split /,/)[0,3,4];
		for ($symbol1,$symbol2){
			$_ = uc($_);
		}

		next if exists $block_locus{$block}{$block_locus_core};
		#next if $block_full_locus_name =~ '^A\d+OA';
		$block_locus{$block}{$block_locus_core} = $symbol1.'-'.$symbol2;
    }
	return(\%block_locus);
}


sub get_WGD_block_Bowers_2003{
	my ($infile, $filter_S) = @_;
	open(my $IN, '<' , "$infile") or die "WGD_block_file cannot be opened";
	my (%strand, %block_locus);
	
	while(<$IN>){
		chomp;
		my (%symbol, $strand, $block_locus_core, $block_locus_num, $block, $symbol);
		# i.e. $_ : >B01fA01#005#000F%1.0007-	At1g01070	>AT1.0007-	A01N006a

		my @line=split /\t/;
		$symbol = uc($line[1]);
		do {last if $line[3] =~ /^S\d+/} if $filter_S;
		#$line[3] =~ /^(\w+)(\w)/;
		($block_locus_core, $block_locus_num) = $line[3] =~ /^(\w+)(\w)/;
		($block) = $block_locus_core =~ /^(.\d\d)/;
		$block_locus{$block}{$block_locus_core} .= $symbol . do{$block_locus_num eq 'a' ? '-' : ''};
	}
	return(\%block_locus);
}

1;
__END__

