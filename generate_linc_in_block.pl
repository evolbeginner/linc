#! /usr/bin/perl

# identify lincs generated via WGD
# Author: Sishuo Wang from Department of Botany, The University of British Columbia
# E-mail: sishuowang@hotmail.ca
# Your help and suggestion are highly appreciated.

# last updated on 2014-12-03
# please use the argument '-h' or '--help' to see the usage

# Note line:
# * next if $seq_name !~ /At\d+NC/;


##############################################################
BEGIN{
    use File::Basename;
	my $dir = dirname($0);
	push @INC, "$dir/module";
}


##############################################################
use 5.010;
use List::Util qw(sum max min first);
use extract_gene;
use Getopt::Long;
use File::Basename;

use get_WGD_block;
no autovivification;


my $extend_gene_width = 10;
my $count_cutoff=10;

my @necessary_params = qw(genome_gff_files linc_gff_files WGD_block_file linc_blast_output block_pair_type);
my (@genome_gff_files, @linc_gff_files, $WGD_block_file, $linc_blast_output, $block_pair_type, $regular_type);
my ($blast_chr_bias_href, $blast_chr_Q_bias, $blast_chr_S_bias, $outdir, $force, $linc_regexp, $is_output_flanking_genes, $no_output_linc, $help);

&show_help() if not @ARGV;


### ---------------------------------------------- ###
$numOfGenesInBlockCutoff = 0;

$regular_type = 'Bowers';
#$linc_blast_output = "../blast/liu_linc_blast_output";
$linc_blast_output = undef;
$blast_chr_Q_bias=undef;
$blast_chr_S_bias=undef;

$linc_blast_e_value_threshold = 1e-5;
$linc_blast_identity_threshold= 80;
$linc_blast_aln_length_threshold= 40;

$linc_regexp = 'At\d+NC|lincRNA|XLOC_|lncRNA';
### ---------------------------------------------- ###


GetOptions(
	'genome_gff=s'		    =>	\@genome_gff_files,
	'linc_gff=s'		    =>	\@linc_gff_files,
	'WGD_block_file=s'	    =>	\$WGD_block_file,
	'linc_blast_output=s'	=>	\$linc_blast_output,
	'block_pair_type=s'	    =>	\$block_pair_type,
    'regular_type=s'        =>  \$regular_type,
    'numOfGenesInBlockMin=s'=>  \$numOfGenesInBlockCutoff,
	'extend_gene_width=s'	=>	\$extend_gene_width,
	'count_cutoff=s'	    =>	\$count_cutoff,
	'A01=s'			        =>	\$A01_fast,
    'e_value|evalue=s'      =>  \$linc_blast_e_value_threshold,
    'identity=s'            =>  \$linc_blast_identity_threshold,
    'length|aln_length=s'   =>  \$linc_blast_aln_length_threshold,
    'blast_chr_Q_bias=s'    =>  \$blast_chr_Q_bias,
    'blast_chr_S_bias=s'    =>  \$blast_chr_S_bias,
    'outdir=s'              =>  \$outdir,
    'force!'                =>  \$force,
    'linc_regexp|lincRNA_regexp=s'  =>  \$linc_regexp,
    'output_flanking'       =>  \$is_output_flanking_genes,
    'no_output_linc!'       =>  \$no_output_linc,
    'h|help!'               =>  \$help,
) or die "illegal params";

&show_help() if $help;

foreach (@necessary_params){
	print "$_ is not defined\n" and &show_help() if not \${$_} and not \@{$_};
	given ($_){
		when (${$_}){
			print $_."\t".${$_}."\n";
		}
		when (@{$_}){
			print $_."\t"."@{$_}"."\n";
		}
	}
}

#$linc_seq_file = "/home/sswang/project/linc/seq/liu_linc.fa";
$block_detailed_step = 50;
#$A01_fast='A15';
undef $filter_S;
undef $adhore;

$blast_chr_bias_href->{'Q'} = $blast_chr_Q_bias if $blast_chr_Q_bias;
$blast_chr_bias_href->{'S'} = $blast_chr_S_bias if $blast_chr_S_bias;

die "linc_blast_output needs to be given! Exiting ......" if not $linc_blast_output;

if (-e $outdir){
    if ($force){
        `rm -rf $outdir` and die "Removal of outdir $outdir failed!";
    }
    else{
        die "Outdir $outdir has already existed!";
    }
}
`mkdir -p $outdir`;


open(my $OUT_params, '>', join('/', $outdir, "params")) || die;
open(my $OUT_printing_pairs, '>', join('/', $outdir, "pairs")) || die;
print $OUT_params "block_detailed_step:\t$block_detailed_step\n";
print $OUT_params "e_value:\t$linc_blast_e_value_threshold\n";
print $OUT_params "identity:\t$linc_blast_identity_threshold\n";
print $OUT_params "aln_length:\t$linc_blast_aln_length_threshold\n";
print $OUT_params "extend_gene_width:\t$extend_gene_width\n";
print $OUT_params "count_cutoff:\t$count_cutoff\n";
print $OUT_params "numOfGenesInBlockMin:\t$numOfGenesInBlockMin\n";
close $OUT_params;

open (my $OUT_linc, ">", join('/', $outdir, "lincs")) || die "" if (! $no_output_linc);
open (my $OUT_all_items, ">", join('/', $outdir, "all_items")) || die "";
open (my $OUT_WGD_pairs, ">", join('/', $outdir, "WGD_pairs")) || die "";


#############################################################################
#############################################################################
while(my $genome_gff_file = shift @genome_gff_files){
	(*posi, *gene_info) = &read_genome_gff_file($genome_gff_file, 'gene', \%gene_info);
}


while (my $linc_gff_file = shift @linc_gff_files){
    *linc_info = &read_linc_gff_file($linc_gff_file);
}
#*linc_info = &extract_gene::read_gene_seq($linc_seq_file, \%linc_info , 'fasta', 'DNA');


(*linc_pair, *linc_paralog_rela) = &read_linc_blast_output($linc_blast_output, 
    $linc_blast_e_value_threshold, $linc_blast_identity_threshold, $linc_blast_aln_length_threshold);

*block = &get_WGD_block($block_pair_type, $WGD_block_file, $filter_S, $regular_type);

*paired_WGD_block_rela = &get_paired_WGD_block(\%block);
print "num of relationships\t" . scalar (keys %paired_WGD_block_rela) . "\n";
output_WGD_pairs(\%paired_WGD_block_rela, $OUT_WGD_pairs);


#### disabled pass ####
*block_info = &generate_block_range(\%block, \%gene_info);
#&print_block_range(\%block_info);

*block_detailed = &generate_detailed_block(\%block_info, $block_detailed_step, \%gene_info);
print "\n";

foreach (qw(gene linc)){
	$hash_name = 'gene_info' if $_ eq 'gene';
	$hash_name = 'linc_info' if $_ eq 'linc';
	(*block_info, *locus_in_block) = 
        &get_block_content(\%{$hash_name}, \%block_detailed, \%block_info, $block_detailed_step, $_, \%locus_in_block);
}

&print_num_within_block(\%locus_in_block);


###############   ********************************************************   ###############
print "\n*************************\n\n";
print "Generating block content info..........................\n";
foreach my $block_name (sort keys %block_info){
	my (%paired_seq_within_block,%seq_order);
    my @sorted_all_seqs;
	my (%a_order);
	my @paired_seq_within_block;
	my @sorted_block_paired;
	my (@num_gene_within_block, @num_linc_within_block);

	if ($A01_fast) {
        next if $block_name !~ /^$A01_fast/;
    }

	print $block_name."\n";
	my $block_name_arrayref = $block_info{$block_name};
	for my $block_num (1..$#{$block_name_arrayref}){ # block_num =1 or 2
		my ($block_name_hashref2, $block_name_hashref3_1, $block_name_hashref3_2);
		my ($block_range, $block_range_pure_num, $block_pos1, $block_pos2, $block_size);
		my (@seq_name);
		my (%block_combined);
		my (%seq_within_block);
		
		print $block_num."\n";
		$block_name_hashref2 = $block_name_arrayref->[$block_num];
		foreach my $locus_type (keys %{$block_name_hashref2}){ # $locus_type = gene or linc
			$block_name_hashref3 = $block_name_hashref2->{$locus_type};
			foreach my $seq_name (keys %{$block_name_hashref3}){
				my ($hash_name, $pos1);
				$hash_name = $locus_type.'_info';
				$pos1 = ${$hash_name}{$seq_name}{pos1}."\n";
				$seq_within_block{$seq_name}{pos1} = $pos1;
				$seq_within_block{$seq_name}{locus_type} = $locus_type;
			}
		}

        @{$sorted_all_seqs[$block_num]}=sort{$seqs{$a}{pos1} <=> $seqs{$b}{pos1}} (keys %seq_within_block);
        # for (@{$sorted_all_seqs[$block_num]}){print $_."\n";}

		foreach (keys %seq_within_block){
			# paired_seq_within_block contais both WGD pairs and lincs
			if ($seq_within_block{$_}{locus_type} eq 'linc'){
				$paired_seq_within_block[$block_num]{$_} = 1;
			}
			if ($seq_within_block{$_}{locus_type} eq 'gene'){
				next if not exists $paired_WGD_block_rela{$_};
				$paired_seq_within_block[$block_num]{$_} = 1;
			}
		}

		my %hash = %{$paired_seq_within_block[$block_num]};
		@{$sorted_block_paired[$block_num]} = 
			sort {$seq_within_block{$a}{pos1} <=> $seq_within_block{$b}{pos1}} keys %hash;
			
        print $OUT_all_items "***\t$block_name\t$block_num\n";
        foreach (@{$sorted_block_paired[$block_num]}){
            print $OUT_all_items $_."\n";
            #print "===\n" if $block_num == 2;
        }

		# only lincs and genes with WGD pairs will be put in @{$sorted_block_paired[$block_num]}

		my @array = @{$sorted_block_paired[$block_num]};
		for (0..$#array){
			my $ele = $array[$_];
			$a_order[$block_num]{$ele} = $_ if $block_num == 1;
			$a_order[$block_num]{$ele} = $_ if $block_num == 2;
			$seq_order[$block_num]{$ele} = $_;
		}
		
		$block_name_hashref3_1 = $block_name_hashref2->{gene};
		$block_name_hashref3_2 = $block_name_hashref2->{linc};
		
		$block_range = $block_name_hashref2->{range}; 		# e.g. 2:1111-9999
		($block_range_pure_num) = $block_range =~ /\:(.+)/; 	#e.g. 1111-9999
		#($block_range_pure_num) = $block_range =~ /^\d+\:(.+)/; 	#e.g. 1111-9999
		($block_pos1,$block_pos2) = split /\-/,$block_range_pure_num;
		$block_size = abs($block_pos1-$block_pos2+1);
		$num_gene_within_block[$block_num] = scalar (keys %$block_name_hashref3_1);
		$num_linc_within_block[$block_num] = scalar (keys %$block_name_hashref3_2);
		#$length_linc_within_block[$block_num] = mean_property($block_name_hashref3_2, \%linc_info, 'length');

		print "block range is\t$block_range\t";
		print "block size is\t$block_size\t";
		print "\n";
		print "Num of genes within block is\t$num_gene_within_block[$block_num]";
		print "\t";
		print "Num of lincs within block is\t$num_linc_within_block[$block_num]";
		print "\n";

        if (! $no_output_linc){
            foreach my $key (keys %$block_name_hashref3_2){
                print $OUT_linc $block_name."-".$block_num."\t".$key."\n";
            }
            print $OUT_linc "\n";
        }
		#print "The average length of lincRNAs in the block is $length_linc_within_block[$block_num]"."\n";
	}

    my ($arrayref1, $arrayref2);
    my $block_num = 1;
    ($arrayref1, $arrayref2) = @sorted_block_paired[1,2];

    &output_adhore("adhore_pair", $arrayref1, $arrayref2) if defined ($adhore);

    foreach my $order1 (0..$#$arrayref1){
        my (@gene1, $gene1_arrayref, $indexs_gene1_aref);
        my ($count, $seq_name, $extend_gene_num_of_max);
        $seq_name = $arrayref1->[$order1];
        $extend_gene_num_of_max = int(2*$extend_gene_width) if not $extend_gene_num_of_max;
        ($gene1_arrayref, $indexs_gene1_aref) = &extend_gene($order1-$extend_gene_width, $order1+$extend_gene_width, $extend_gene_num_of_max, 0, $#$arrayref1, $arrayref1, 'array');
        @gene1 = @$gene1_arrayref;

        # -------------------------------------------------- #
        #next if $seq_name !~ /$linc_regexp/i;
        next if not exists $linc_info{$seq_name};
        # -------------------------------------------------- #

        if (exists $linc_paralog_rela{$seq_name}){
            #map {print $_."\n"} keys %{$linc_paralog_rela{$seq_name}};
            my $element_hash = $linc_paralog_rela{$seq_name};
            for my $paralog (keys %$element_hash){
                my @paired_genes;
                my ($count,$count2,%gene2);
                $count = 0;
                next if not exists $paired_seq_within_block[3-$block_num]{$paralog};
                my $order2 = $a_order[3-$block_num]{$paralog}+1;
                my ($gene2_hashref, $indexs_gene2_aref) = &extend_gene( $order2-$extend_gene_width, 
                    $order2+$extend_gene_width, $extend_gene_num_of_max,
                    0, $#$arrayref2,
                    $arrayref2, 'hash' );
                my %gene2 = %$gene2_hashref;
                my $count_cutoff = int((scalar(keys %gene2) + scalar(@gene1))/4);
                for my $gene (@gene1){
                    foreach my $paralog (keys %{$paired_WGD_block_rela{$gene}}){
                        my $WGD_paired_gene = $paralog;
                        if (exists $gene2{$WGD_paired_gene}){
                            ++$count;
                            push @paired_genes, join("\t", $gene, $WGD_paired_gene);
                            last;
                        }
                    }
                }
                next if $count < $count_cutoff;
                push @paired_genes, join("\t",$seq_name,$paralog);

                # ********************************************************************
                if ($is_output_flanking_genes){
                    map {print $_."\n"} @paired_genes;
                    map {print $_."\n"} @gene1;
                    map {print $_."\n"} keys %gene2;
                }
                # ********************************************************************

                #print "$indexs_gene1_aref->[0]\t$indexs_gene1_aref->[1]\n";
                #my $index_aref = &get_index($sorted_all_seqs[$block_num],[@$arrayref1[@$indexs_gene1_aref]]);
                #print join("\t",@$index_aref)."\n";
                #for my $index ($indexs_gene2_aref->[0]..$indexs_gene2_aref->[1]){print $gene2{$index}."\n";}

                my $pair_name = join ("-", sort $seq_name,$paralog);
                my $e_value = $linc_pair{$pair_name}{e_value};
                my $printing_item = join ("\t",$block_name,$seq_name,$paralog,$count,
                    scalar @gene1, scalar keys %gene2,
                    $e_value,@num_gene_within_block );
                print $printing_item."\n";
                print $OUT_printing_pairs $printing_item."\n";
            }
        }
    }
}

close  $OUT_printing_pairs;
close $OUT_linc if (! $no_output_linc);



##################################################################################
sub extend_gene{
	my ($start, $end, $num_of_max, $lower, $upper, $arrayref1, $type) = @_;
	my $return_ref;
    my $middle_point = ($start+$end)/2;
    my $counter_limit = ($end-$start)/2;
    my ($counter,$counter2)=(1,0);

    while($middle_point-$counter>=$lower){
        last if $counter2 >= $counter_limit;
        ++$counter;
        $a=$middle_point-$counter;
        #next if $arrayref1->[$a] =~ /$linc_regexp/i;
        next if exists $linc_info{$arrayref1->[$a]};
        my $neighbor_seq_name = $arrayref1->[$a];
        if ($type eq 'hash'){
                $return_ref->{$neighbor_seq_name} = 1 if scalar keys %$return_hashref<$num_of_max;}
        elsif ($type eq 'array'){
                push @$return_ref, $neighbor_seq_name if scalar @$return_ref<$num_of_max;
        }
        $counter2++;
    }

    my $lower_counter=$middle_point-$counter+1;
    my ($counter,$counter2)=(1,0);
    #print join("\t",$start,$end,$counter,$counter_limit,$middle_point)."\n";

    while($middle_point+$counter<=$upper){
        last if $counter2 >= $counter_limit;
        ++$counter;
        $a=$middle_point+$counter;
        #next if $arrayref1->[$a] =~ /$linc_regexp/i;
        next if exists $linc_info{$arrayref1->[$a]};
        my $neighbor_seq_name = $arrayref1->[$a];
        if ($type eq 'hash'){
                $return_ref->{$neighbor_seq_name} = 1 if scalar keys %$return_hashref<$num_of_max;
        } elsif ($type eq 'array'){
                push @$return_ref, $neighbor_seq_name if scalar @$return_ref<$num_of_max;
        }
        $counter2++;
    }
    my $upper_counter=$middle_point+$counter-1;
	return ($return_ref, [$lower_counter,$upper_counter]);
}



sub get_WGD_block{
	my ($block_pair_type, $WGD_block_file, $filter_S, $regular_type) = @_;
	given ($block_pair_type){
		when($_ eq 'Bowers')	{(*block) = &get_WGD_block_Bowers_2003($WGD_block_file, $filter_S);}
		when($_ eq 'regular')	{(*block) = &get_WGD_block_regular($WGD_block_file, $filter_S, $regular_type, $numOfGenesInBlockCutoff);}
		default                 {die "Illegal para for block_pair_type has been given"}
	}
	return (\%block);
}


sub get_block_content
{
    my ($gene_info_ref, $block_detailed_ref, $block_info_ref, $block_detailed_step, $type, $locus_in_block_ref) = @_;

    print "Getting block content for $type ......\n";
    # type = 'gene' or 'linc'
    my %gene_info = %$gene_info_ref;
    my %block_detailed = %$block_detailed_ref;
    my $block_info = %$block_info_ref;
    my $locus_in_block = %$locus_in_block_ref;
    my $total = scalar keys (%gene_info);
    my $counter = 0;
    my $mei_fen = int($total/20);

    foreach my $gene_name (keys %gene_info){
        print "." if ++$counter % $mei_fen == 0;
        my ($posi, $chr, $strat, $end, @start_array);
        $posi = $gene_info{$gene_name}{posi};
        # ($chr) = $gene_name =~ /^\D*(\d+)/;
        $chr = $gene_info{$gene_name}{'chr'};
        ($start, $end) = split /\-/, $posi;
        do {push @start_array, $start+$_;} for (0..$block_detailed_step-1);

        while (@start_array){
            my $new_start = shift @start_array;
            my $start_hashref = $block_detailed{$chr};
            if (exists $start_hashref -> {$new_start}){
                $linc_SuoShu_block_hashref = $start_hashref -> {$new_start};
                for my $linc_SuoShu_block (keys %{$linc_SuoShu_block_hashref}){
                    #linc_SuoShu_block refers to the block that a lincRNA belongs to
                    #my $linc_SuoShu_block = $array_ref->{$new_start};
                    $locus_in_block{$type}{$gene_name} = "";
                    ($block_name, $block_num) = $linc_SuoShu_block =~ /^(.+)\:(\d+)$/;
                    $block_info{$block_name}[$block_num]{$type}{$gene_name} = "";
                }
                last;
            }
        }
    }
    print "\n";
    return (\%block_info, \%locus_in_block);
}


sub generate_block_range
{
    my ($block_ref, $gene_info_href) = @_;
    my %block = %$block_ref;
    my $block_printing_counter; # print_block_counter

    my $num_of_keys_of_hash_block = scalar keys %block;
    print "Generating block range\n";

    my $block_printing_countersh_href = &generate_printing_countersh($num_of_keys_of_hash_block);

    iterate_block: foreach my $block_name (sort {$a<=>$b} keys %block){
        if ($A01_fast){
            next if $block_name !~ /^$A01_fast/;
        }
        $block_printing_counter++ ;
        system "echo -ne '.'" if exists $block_printing_countersh_href->{$block_printing_counter};
        my (@coordinate, %array2, @chr_num, @chr_num_tmp);
        my $hash2_ref = $block{$block_name};
        foreach $key2 (sort keys %$hash2_ref){
            my $res = $hash2_ref->{$key2};
            ($gene1,$gene2) = split /[!]/,$res;

            # e.g., Chr1 or chr1 -> 1
            $chr_num_tmp[1] = $gene_info_href->{$gene1}{'chr'} if exists $gene_info_href->{$gene1};
            $chr_num_tmp[2] = $gene_info_href->{$gene2}{'chr'} if exists $gene_info_href->{$gene2};
            #($chr_num[1]) = $gene1 =~ /(\d+)/;
            #($chr_num[2]) = $gene2 =~ /(\d+)/;

            foreach (1..2){
                if (not $chr_num[$_]){
                    $chr_num[$_] = $chr_num_tmp[$_];
                }
                else{
                    if ($chr_num[$_] ne $chr_num_tmp[$_]){
                        print "chr_num inconsistent in $block_name\t$_\n";
                        next iterate_block;
                    }
                }
            }

            foreach my $key3(1..2){
                foreach my $key4(1..2){
                    #$coordinate[$key3][$key4] .= $posi{${'gene'.$key3}}{'pos'.$key4}." ";
                    $coordinate[$key3][$key4] .= $gene_info{${'gene'.$key3}}{'pos'.$key4}." ";
                }
            }
        }

        foreach my $key1(1..2){
            my (%array2,%max,%min);
            foreach my $key2(1..2){
                my @a = split /\s+/, $coordinate[$key1][$key2];
                $max[$key2] = max @a;
                $min[$key2] = min @a;
            }

            foreach my $key2 (1..2){ 
                my $max=$max[$key2];
                my $min=$min[$key2];
                $array2{$max.'-'.$min} = abs($max-$min+1);
            }
            
            my $final_max = max values %array2;
            $final_posi = first {$array2{$_} == $final_max} keys %array2; 
            $block_info{$block_name}[$key1]{range} = "$chr_num[$key1]:$final_posi";
            #print $block_info{$block_name}[$key1]{range}."\n";
        }
    }
    print "\n";
    return(\%block_info);
}


sub generate_detailed_block
{
    my ($block_info_ref, $step, $gene_info_hashref) = @_;
    my (%range, $block_info, %block_detailed, $num_of_keys_of_hash_block_info);
    my %block_info = %$block_info_ref;

    $num_of_keys_of_hash_block_info = scalar keys %block_info;
    %gene_info = %$gene_info_hashref;
    my $block_printing_countersh_href = &generate_printing_countersh($num_of_keys_of_hash_block_info);

    print "\nGenerating detailed block\n";
    #&print_dot($num_of_keys_of_hash_block_info);
    for my $block_name(sort keys %block_info){
        do {next if $block_name ne $A01_fast} if $A01_fast;
        my (@range,$hash2_ref,$range,$chr);
        $block_printing_counter++ ;
        system "echo -ne '.'" if exists $block_printing_countersh_href->{$block_printing_counter};
        $hash2_ref = $block_info{$block_name};
        for my $key1 (1..2){
            #($range = $block_info{$block_name}[$key1]{range}) =~ /^\D*(\d+)\:([^\-]+)\-([^\-]+)/;
            ($range = $block_info{$block_name}[$key1]{range}) =~ /^([^:]+)\:([^\-]+)\-([^\-]+)/;
            $chr = $1;
            @range = do {$2<$3 ? ($2,$3) : ($3,$2)};
            for my $coordinate ($range[0]..$range[1]){
                if ($coordinate % $step == 0){
                    #$block_detailed[$chr]{$coordinate} = $block_name.":$key1";
                    $block_detailed{$chr}{$coordinate}{$block_name.":$key1"}=1;
                }
            }
        }
    }
    print "\n";
    return (\%block_detailed);	
}


sub get_paired_WGD_block{
	my ($block_ref) = @_;
	my %paired_WGD_block_rela;
	for my $block_name (keys %$block_ref){
		my $hashref1 = $block_ref->{$block_name};
		get_paired_WGD_block__cycle1: for my $block_locus_name (keys %$hashref1){
			my $symbol_combined = $hashref1->{$block_locus_name};
			my ($symbol1, $symbol2) = split /\!/, $symbol_combined;
			#for (($symbol1,$symbol2)){next get_paired_WGD_block__cycle1 if exists $paired_WGD_block_rela{$_};}
			#$paired_WGD_block_rela{$symbol1} = $symbol2;
			#$paired_WGD_block_rela{$symbol2} = $symbol1;
			$paired_WGD_block_rela{$symbol1}{$symbol2} = '';
		    $paired_WGD_block_rela{$symbol2}{$symbol1} = '';
		}
	}
	return (\%paired_WGD_block_rela);
}


sub print_block_range{
	my ($block_info_ref) = @_;
	print "\nPrinting block_range...............\n";
	my %block_info = %{$block_info_ref};
	foreach my $key1 (sort keys %block_info){
		for (1..2){
                    print $key1."\t".$block_info{$key1}[$_]{range}."\n"
                }
	}
}


sub print_num_within_block{
    my ($locus_in_block_ref) = @_;
    my %locus_in_block = %$locus_in_block_ref;
    print "\nThe number of genes and lincRNAs that are localized in WGD blocks are:\n";
    foreach my $key1 (keys %locus_in_block){
        print $key1. "\t";
        print scalar keys %{$locus_in_block{$key1}};
        print "\n";
    }
}


sub output_WGD_pairs{
    my ($paired_WGD_block_rela_href, $OUT_WGD_pairs) = @_;
    foreach my $key (keys %$paired_WGD_block_rela_href){
        my @a = keys(%{$paired_WGD_block_rela{$key}});
        print $OUT_WGD_pairs $key."\t".$a[0]."\n";
    }
    close $OUT_WGD_pairs;
}


##################################################################################
sub read_linc_gff_file{
    # strand is not considered here
	my ($gff_file) = @_;
	my %linc_info;
	open (my $IN, '<', $gff_file) or die "linc gff file cannot be opened";
	while(<$IN>){
		chomp;
		my ($chr,$name,$start,$end) = split /\t/;
        $linc_info{$name}{'chr'} = $chr;
		#my ($chr_num) = $chr =~ /(\d+)/;
        #$linc_info{$name}{'chr'} = $chr_cum;
		$linc_info{$name}{posi} = $start.'-'.$end;
		$linc_info{$name}{pos1} = $start; # pos1 is not posi 
	}
	close $IN;
	return (\%linc_info);
}


sub read_genome_gff_file{
	my ($gff_file, $type0, $input_gene_info_hashref) = @_;
	my (%posi, %gene_info);
	%gene_info = %$input_gene_info_hashref if $input_gene_info_hashref;

	open ($IN, '<' , "$gff_file") or die "gff_file $gff_file cannot be opened";
	while(<$IN>){
		chomp;
		my ($chr, $type, $pos1, $pos2, $strand, $length0, $length, $ref2);
		($chr, $type, $pos1, $pos2, $strand, $attributes) = (split /\t/)[0,2,3,4,6,8];
        #($chr) = $chr =~ /(\d+)/; # e.g., chr5 -> 5
		next if $type ne $type0;
		
		#($gene) = $attributes =~ /^ID\=([.A-Za-z0-9]+)/;
		if ($attributes =~ /ID\=([^;]+)/){
            $ref2 = $posi{$1};
            $gene = $1;
        }
        elsif ($attributes =~ /Name\=([^;]+)/){
            $ref2 = $posi{$1};
            $gene = $1;
        }
        $gene = uc($gene);

		$length = abs($pos1-$pos2+1);
		if (defined $ref2){
			$length0 = abs($ref2->{pos1} - $ref2->{pos2} + 1);
			next if $length0 > $length;
		}
		$gene_info{$gene}{posi} = $pos1.'-'.$pos2;
		$gene_info{$gene}{pos1} = $pos1;
		$gene_info{$gene}{pos2} = $pos2;
		$gene_info{$gene}{chr} = $chr;
	}	
	return(\%posi, \%gene_info);
}


sub read_linc_blast_output
{
	my ($blast_output, $e_value_threshold, $identity_threshold, $aln_length_threshold) = @_;
	my (%pair,%para_rela);
	$e_value_threshold = 1e-5 if not defined $e_value_threshold;
    $identity_threshold= 80 if not defined $identity_threshold;
    $aln_length_threshold= 40 if not defined $aln_length_threshold;

	open (my $IN, '<', $blast_output) or warn "blast_output file cannot be opened";
	while(<$IN>){
		chomp;
		my ($seq1, $seq2, $identity, $aln_length, $e_value) = (split /\t/)[0,1,2,3,10];
		next if $seq1 eq $seq2;
        ### ------------------------------------------------- ###
        if ($e_value > 1e-10){
		    if ($e_value <= $e_value_threshold){
                next if $identity < $identity_threshold or $aln_length < $aln_length_threshold;
            }
            else{
                next;
            }
        }
		my $pair_name = join ("-", sort ($seq1, $seq2));
		$pair{$pair_name}{e_value} = $e_value if (not exists $pair{$pair_name} or $pair{$pair_name}{e_value} > $e_value);
		# $pair{$pair_name}{identity} = $identity if (not exists $pair{$pair_name} or $pair{$pair_name}{e_value} > $e_value);
		$para_rela{$seq1}{$seq2} = 1;
		$para_rela{$seq2}{$seq1} = 1;
	}
	close $IN;
	return (\%pair,\%para_rela);
}


sub get_index{
    my ($aref, $ele_aref) = @_;
    my %hash;
    my @indexs;
    my @eles=@$ele_aref;
    map {$hash{$aref->[$_]}=$_} 0..$#$aref;
    for my $ele (@eles){
        push @indexs, $hash{$ele};
    }
    return(\@indexs);
}


########################################################################
sub print_dot{
	my ($num_of_dots) = @_;
	print '.' x $num_of_dots;
	print "\n";
}


sub generate_printing_countersh{
    my %block_printing_countersh;
    my ($num_of_keys_of_hash_block) = @_;

    if ($num_of_keys_of_hash_block <= 20){
        map { $block_printing_countersh{$_}='' } 1..$num_of_keys_of_hash_block;
    } else{
        my $num_per_fen = $num_of_keys_of_hash_block/20;
        foreach (1..20){
            my $num = int($_*$num_per_fen);
            $block_printing_countersh{$num} = '';
        }
    }

    return (\%block_printing_countersh);
}


sub show_help{
    print "perl " . basename($0) . "\n";
    print <<EOF
[Required arguments]
--genome_gff
--linc_gff
--WGD_block_file
--linc_blast_output

[Optional arguments]
--block_pair_type
--extend_gene_width
--count_cutoff
--A01
--outdir
--e_value|evalue
--identity
--length|aln_length
--force             default: off
--output_flanking
--no_output_linc    default: off
--h|help

EOF
;
    print 'Please write e-mails to sishuowang@hotmail.ca if you have any question and/or suggestion. Your help is highy appreciated!' . "\n";
    print "\n";
    exit();
}


sub mean{
	return 0 if scalar @_ == 0;
	return (sum (@_)/@_);
}



#############################################################################
# finalize associating
sub mean_property{
	my ($hashref3_2, $info_hashref, $property) = @_;
	my (@property, $mean_length);
	foreach (keys %$hashref3_2){
		push @property, $info_hashref->{$_}{$property};
	}
	$mean_length = mean (@property);
	return ($mean_length);
}


sub output_adhore
{
	my ($outfile, $arrayref1, $arrayref2) = @_;	
	my $kkk;
	open (my $pair_OUT, '>', "$outfile");
	foreach my $name_ref (($arrayref1,$arrayref2)){
		$kkk++;
		open (my $OUT, '>', "pair_rela_quan.$kkk");
		foreach my $seq_name (@$name_ref){
			print $OUT "$seq_name+\n";
			if (exists $linc_paralog_rela{$seq_name}){
				for (keys %{$linc_paralog_rela{$seq_name}}){
					print $pair_OUT "$seq_name\t$_\n";
				}
			}
		}
		close $OUT;
	}
	
	foreach my $symbols (values %{$block{$block_name}}){
		my $out;
		($out = $symbols) =~ s/\-/\t/;
		print $pair_OUT $out."\n";
	}
	close $pair_OUT;
}




################################################################################
#          *************************** ·ÏÎï **************************         #
################################################################################
sub get_WGD_block_Bowers_2004{
	my ($infile, $filter_S) = @_;
	open($IN, '<' , "$infile") or die "WGD_block_file cannot be opened";
	my (%strand, %block_locus);
	
	while(<$IN>){
		chomp;
		my (%symbol, $strand, $block_locus_core, $block_locus_num, $block);
		# e.g. $_ : >B01fA01#005#000F%1.0007-	At1g01070	>AT1.0007-	A01N006a

		@line=split /\t/;
		$symbol = uc($line[1]);
		#($strand) = $line[2] =~ /([-+])$/; $strand{$symbol} = $strand;
		do {last if $line[3] =~ /^S/} if defined $filter_S;
		($block_locus_core, $block_locus_num) = $line[3] =~ /^(\w+)(\w)/;
		($block) = $block_locus_core =~ /^(.\d\d)/;
		$block_locus{$block}{$block_locus_core} .= $symbol . do {$block_locus_num eq 'a' ? '-' : ''};
	}
	return(\%block_locus);
}


sub USAGE{
	print "USAGE:\tperl $0 ";
	print "<genome_gff=\$genome_gff_file> <linc_gff=\$lin_gff_file> <WGD_block_file=\$WGD_block_file> [linc_blast_output=\$linc_blast_output] \n\n";
	exit();
}


#############################################################################
=cut
$genome_gff_file = '../data/TAIR/TAIR10/TAIR10_GFF3_genes.gff';
$linc_gff_file   = '../data/pseudo/w7/liu_linc.ATH.blast.out.gff.pseudo';
$WGD_block_file  = '../data/Bowers_block.list.2003';
$linc_blast_output = '../data/pseudo/w7/liu_linc.ATH.blast.out.pseudo';
$block_pair_type = 'Bowers';
=cut


