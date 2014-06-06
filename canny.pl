#!/usr/bin/perl

use strict;
use warnings;
use Vcf;
use Statistics::R;
use Getopt::Long;

my $version = "0.0.6";

#Version
#   0.0.6
#       -added option to input range as alternative to filename
#
#   0.0.5
#       -added option to output gentypes in CN format
#       -added upper bound on cgh reference sequence correction
#
#   0.0.4
#       -added option to use CGH reference sequence to correct cluster assignment
#
#   0.0.3
#       -added option to filter sites with largest cluster exceeding std limits
#       -added option to use previously calculated cluster directory
#       -changed plotting options and format
#
#   0.0.2
#       -Modified cluster assignment to genotypes
#       -Added Getopt module for parameter handling

my %opts = ();
$opts{min_mapq}        = 10;
$opts{read_len}        = 101;
$opts{bandwidth}       = 0.05;
$opts{threshold}       = 0.3;
$opts{min_num_probes}  = 5;
$opts{max_cluster_std} = 0.5;
$opts{verbose}         = 0;
$opts{work_dir}        = ".";
$opts{ratio_filename}  = "/scratch/remills_flux/remills/1kg_aCGH/ratios/allData.filtered.gc.txt.corr3.allCols.dat.gz";
$opts{sample_filename} = "/scratch/remills_flux/remills/1kg_aCGH/samples.txt";
$opts{cluster_label}   = "cluster";
$opts{samtools}        = "samtools";
$opts{output_format}   = "CN";
$opts{version} = $version;

my $optResult = GetOptions(
    "input_filename=s"        => \$opts{input_filename},
    "input_range=s"           => \$opts{input_range},
    "output_filename=s"       => \$opts{output_filename},
    "ratio_filename=s"        => \$opts{ratio_filename},
    "sample_filename=s"       => \$opts{sample_filename},
    "work_dir=s"              => \$opts{work_dir},
    "bandwidth=f"             => \$opts{bandwidth},
    "threshold=f"             => \$opts{threshold},
    "cluster_label=s"         => \$opts{cluster_label},
    "min_num_probes=i"        => \$opts{min_num_probes},
    "max_cluster_std=f"       => \$opts{max_cluster_std},
    "plot_results"            => \$opts{plot_results},
    "use_existing_clusters=s" => \$opts{use_existing_clusters},
    "cgh_ref_seq=s"           => \$opts{cgh_ref_seq},
    "cgh_ref_mean_cov=f"      => \$opts{cgh_ref_mean_cov},
    "read_len=i"              => \$opts{read_len},
    "min_mapq=i"              => \$opts{min_mapq},
    "samtools=s"              => \$opts{samtools},
    "output_format=s"         => \$opts{output_format},
    "verbose"                 => \$opts{verbose}
);

checkOptions( $optResult, \%opts, $opts{version} );

if ( !-d "$opts{work_dir}/clusters" && !$opts{use_existing_clusters} ) { mkdir "$opts{work_dir}/clusters"; }
if ( !-d "$opts{work_dir}/results" ) { mkdir "$opts{work_dir}/results"; }
if ( !-d "$opts{work_dir}/medians" ) { mkdir "$opts{work_dir}/medians"; }
if ( !-d "$opts{work_dir}/plots" && $opts{plot_results} ) { mkdir "$opts{work_dir}/plots"; }

my @samples = ();
open( NM, $opts{sample_filename} ) || die "Can't open $opts{sample_filename}, $!\n";
if ( $opts{output_filename} ) {
    open( RS, "| bgzip -c > $opts{work_dir}/results/$opts{output_filename}" ) || die "Can't open $opts{work_dir}/results/$opts{output_filename} for output$!\n";
}
else {
    open( RS, ">&", \*STDOUT ) or die;
}
print RS <<HEADER;
##fileformat=VCFv4.1
##INFO=<ID=LDAF,Number=1,Type=Float,Description="MLE Allele Frequency Accounting for LD">
##INFO=<ID=AVGPOST,Number=1,Type=Float,Description="Average posterior probability from MaCH/Thunder">
##INFO=<ID=RSQ,Number=1,Type=Float,Description="Genotype imputation quality from MaCH/Thunder">
##INFO=<ID=ERATE,Number=1,Type=Float,Description="Per-marker Mutation rate from MaCH/Thunder">
##INFO=<ID=THETA,Number=1,Type=Float,Description="Per-marker Transition rate from MaCH/Thunder">
##INFO=<ID=CIEND,Number=2,Type=Integer,Description="Confidence interval around END for imprecise variants">
##INFO=<ID=CIPOS,Number=2,Type=Integer,Description="Confidence interval around POS for imprecise variants">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description="Length of base pair identical micro-homology at event breakpoints">
##INFO=<ID=HOMSEQ,Number=.,Type=String,Description="Sequence of base pair identical micro-homology at event breakpoints">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=AC,Number=.,Type=Integer,Description="Alternate Allele Count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total Allele Count">
##INFO=<ID=NPROBES,Number=1,Type=Integer,Description="Number of aCGH Probes in region">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=CNV,Description="Copy Number Variation">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="Genotype dosage from MaCH/Thunder">
##FORMAT=<ID=GL,Number=.,Type=Float,Description="Genotype Likelihoods">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele, ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/technical/reference/ancestral_alignments/README">
##INFO=<ID=AF,Number=1,Type=Float,Description="Global Allele Frequency based on AC/AN">
##INFO=<ID=AMR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AMR based on AC/AN">
##INFO=<ID=ASN_AF,Number=1,Type=Float,Description="Allele Frequency for samples from ASN based on AC/AN">
##INFO=<ID=AFR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from AFR based on AC/AN">
##INFO=<ID=EUR_AF,Number=1,Type=Float,Description="Allele Frequency for samples from EUR based on AC/AN">
##INFO=<ID=VT,Number=1,Type=String,Description="indicates what type of variant the line represents">
##INFO=<ID=SNPSOURCE,Number=.,Type=String,Description="indicates if a snp was called when analysing the low coverage or exome alignment data">
##reference=GRCh37
HEADER
print RS "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
while (<NM>) {
    chomp;
    push @samples, $_;
    print RS "\t$_";
}
print RS "\n";
close NM;

my $vcf;
my $x;

if ( $opts{input_filename} ) {
    $vcf = Vcf->new( file => $opts{input_filename} );
    $vcf->parse_header();
}
elsif ( $opts{input_range} ) {
    my ( $chr, $pos, $end ) = $opts{input_range} =~ /(\S+):(\d+)-(\d+)/;
    $vcf = Vcf->new();
    $vcf->add_columns();
    $vcf->add_header_line( { key => 'ALT',    ID => 'CNV', Description => "Query CNV Region" } );
    $vcf->add_header_line( { key => 'INFO',   ID => 'END', Number      => 1, Type => 'Integer', Description => "End coordinate of this variant" } );
    $vcf->add_header_line( { key => 'INFO',   ID => 'NPROBES', Number  => 1, Type => 'Integer', Description => "Number of aCGH Probes in region" } );
    $vcf->add_header_line( { key => 'FORMAT', ID => 'GT',  Number      => 1, Type => 'String', Description => "Genotype" } );
    $vcf->add_header_line( { key => 'FORMAT', ID => 'CN',  Number      => 1, Type => 'String', Description => "Copy number genotype for imprecise events" } );
    $vcf->format_header();

    $$x{CHROM}  = $chr;
    $$x{POS}    = $pos - 1;
    $$x{ID}     = '.';
    $$x{ALT}    = ["<CNV>"];
    $$x{REF}    = 'N';
    $$x{QUAL}   = '.';
    $$x{FILTER} = ['.'];
    $$x{INFO}   = { END => $end };
}
if ( defined( $$x{ID} ) ) {
    clusterRegion($x);
}
else {
    while ( $x = $vcf->next_data_hash() ) {
        clusterRegion($x);
    }
}
close RS;
$vcf->close();

sub clusterRegion {
    my ($x) = @_;
    warn "Analyzing $$x{ID}...\n" if $opts{verbose};
    my $id  = $$x{ID};
    my $chr = $$x{CHROM};
    my $pos = $$x{POS} + 1;
    my $end = $$x{INFO}{END};

    if ( $id eq "." ) { $id = "$chr\_$pos\_$end"; }

    my @ratios = ();
    open( PRB, "tabix $opts{ratio_filename} chr$chr:$pos-$end |" ) || die "$!\n";
    my $probeNum = 0;
    while (<PRB>) {
        my @row = split(/\t/);
        push @ratios, [ splice( @row, 7 ) ];    #column 8 onwards, according to Ankit....
        $probeNum++;
    }

    if ( $probeNum < $opts{min_num_probes} ) { return; }
    my @medians   = ();
    my $sampleNum = 0;
    while ( defined( $ratios[0][$sampleNum] ) ) {
        my @row = ();
        for ( my $i = 0 ; $i < $probeNum ; $i++ ) {
            push @row, $ratios[$i][$sampleNum];
        }
        push @medians, median(@row);
        $sampleNum++;
    }

    open( OUT, ">$opts{work_dir}/medians/$id.dat" );
    print OUT "sample\tmedian\n";
    for ( my $i = 0 ; $i <= $#medians ; $i++ ) {
        print OUT "$samples[$i]\t$medians[$i]\n";
    }
    close OUT;

    my $d;
    my $m;

    if ( $opts{use_existing_clusters} ) {
        open( CL, "$opts{use_existing_clusters}/$id\_clustered.txt" ) || die "Could not open $opts{work_dir}/clusters/$id\_clustered.txt for input with use_existing_clusters, exiting\n";
        while (<CL>) {
            chomp;
            my @row = split(/\t/);
            push @$d, @row;
        }
        close CL;

        open( ML, "$opts{use_existing_clusters}/$id\_means.txt" ) || die "Could not open $opts{work_dir}/clusters/$id\_means.txt for input with use_existing_clusters, exiting\n";
        while (<ML>) {
            chomp;
            my @row = split(/\t/);
            push @$m, @row;
        }
    }
    else {
        my $R = Statistics::R->new();
        $R->startR;
        $R->run(qq`library(LPCM)`);
        $R->run(qq`library(ggplot2)`);
        $R->run(qq`a<-read.table("$opts{work_dir}/medians/$id.dat", header=T, row.names="sample")`);
        $R->run(qq`fit<-ms(a, h=$opts{bandwidth}, scaled=F, thr=$opts{threshold})`);
        $R->run(qq`d<-data.frame(fit\$data)`);
        $R->run(qq`m<-data.frame(fit\$cluster.center)`);
        $R->run(qq`d\$class <- fit\$$opts{cluster_label}.label`);

        if ( $opts{plot_results} ) {
            $R->run(qq`pdf("$opts{work_dir}/plots/$id.pdf")`);

            #$R->run(qq`ggplot(d,aes(x=median, fill=as.factor(class))) + geom_histogram(aes(y=..density..), binwidth=0.01) + geom_density(alpha=0.2) + guides(fill=guide_legend(title="clusters"))`); #include density plot
            $R->run(qq`ggplot(d,aes(x=median, fill=as.factor(class))) + geom_histogram(binwidth=0.01) + guides(fill=guide_legend(title="clusters"))`);    #histogram counts only
            $R->run(qq`dev.off()`);
        }
        $R->run(qq`write.table(d, "$opts{work_dir}/clusters/$id\_clustered.txt", sep="\t", quote=F)`);
        $R->run(qq`write.table(m, "$opts{work_dir}/clusters/$id\_means.txt", sep="\t", quote=F)`);
        $d = $R->get(qq`d`);
        $m = $R->get(qq`m`);
        $R->stopR();
    }

    #if ($#$m > 6) { next; } #TRY
    if ( !clusterQC($d) ) { next; }

    my %key;
    getGenoKey( \%key, $m, $d, $chr, $pos, $end );

    my $genoLine = "";
    my $ac       = 0;
    my $an       = 0;
    for ( my $i = 2 ; $i <= $#$d ; $i++ ) {
        my $sample = $$d[$i];
        $i++;
        my $median = $$d[$i];
        $i++;
        my $index = $$d[$i];
        $genoLine .= "\t$key{$index}{geno}";
        if ( $key{$index}{geno} !~ /\./ ) {
            if ( $opts{output_format} eq "GL" ) {
                my ($cnt) = $key{$index}{geno} =~ tr/1//;
                $ac += $cnt;
                $an += 2;
            }
            elsif ( $opts{output_format} eq "CN" ) {
                $ac += abs( 2 - $key{$index}{geno} );
            }
            $an += 2;
        }
    }

    my $infoLine = format_line_hash( $x, $ac, $an, $probeNum );
    print RS "$infoLine$genoLine\n";
}

sub clusterQC {

    #checks if stdev of largest cluster is within parameter bounds
    my ($d) = @_;

    my $max_count      = 0;
    my $max_mean_index = 0;
    my %medians        = ();
    my %counts         = ();

    #calculate cluster countes
    for ( my $i = 2 ; $i <= $#$d ; $i++ ) {
        my $sample = $$d[$i];
        $i++;
        my $median = $$d[$i];
        $i++;
        my $cluster = $$d[$i];
        $counts{$cluster}++;
        push @{ $medians{$cluster} }, $median;

        if ( $counts{$cluster} > $max_count ) {
            $max_count      = $counts{$cluster};
            $max_mean_index = $cluster;
        }
    }
    my $std = stdev( @{ $medians{$max_mean_index} } );
    if   ( $std <= $opts{max_cluster_std} ) { return 1; }
    else                                    { return 0; }
}

sub getGenoKey {
    my ( $key, $m, $d, $chr, $pos, $end ) = @_;

    my @exp            = qw(-3 -2 -1 0 0.6 1 2);         #based on expectations of log2 ratios per copy number
    my @genos          = qw(1/1 1/1 0/1 0/0 ./. ./.);    #bi-allelic for now, need to update to multi-allelic format for VCF4.1
    my %means          = ();
    my %assigned       = ();
    my %counts         = ();
    my $max_count      = 0;
    my $max_mean_index = 0;

    #calculate cluster countes
    for ( my $i = 2 ; $i <= $#$d ; $i++ ) {
        my $sample = $$d[$i];
        $i++;
        my $median = $$d[$i];
        $i++;
        my $cluster = $$d[$i];
        $counts{$cluster}++;
        if ( $counts{$cluster} > $max_count ) {
            $max_count      = $counts{$cluster};
            $max_mean_index = $cluster;
        }
    }

    #calculate distance from 0
    for ( my $i = 0 ; $i <= $#$m ; $i++ ) {
        if ( $$m[$i] eq 'fit.cluster.center' ) { next; }
        my $index = $$m[$i];
        $i++;
        my $mean = $$m[$i];
        $$key{$index}{mean} = $mean;
    }

    my @sorted = sort { $$key{$a}{mean} <=> $$key{$b}{mean} } keys %$key;

    #find closest cluster to 0
    my $max_index = 0;
    my $min       = abs( $$key{ $sorted[0] }{mean} );
    for ( my $i = 1 ; $i <= $#sorted ; $i++ ) {
        if ( abs( $$key{ $sorted[$i] }{mean} ) < $min ) {
            $min       = abs( $$key{ $sorted[$i] }{mean} );
            $max_index = $i;
        }
    }

    #set largest cluster to 0 as 0/0, and other clusters incremental from there
    my $index = 0;
    for ( my $i = 0 ; $i <= $#sorted ; $i++ ) {
        if ( $sorted[$i] == $max_mean_index ) {
            $index = $i;
            last;
        }
    }

    my $zero = 3;
    my $cn   = 2;
    my $skip = 0;
    if ( $opts{cgh_ref_seq} ) {
        #####update "zeroed" index by taking CGH reference sequence coverage into account - not yet, will try later
        #skip if reference sequence coverage is signifcantly lower than expected - for now
        my $mean_cov = getMeanCov( $chr, $pos, $end );
        my $pObs = 0;
        for ( my $i = 0 ; $i <= $mean_cov ; $i++ ) {
            $pObs += poisson( $i, $opts{cgh_ref_mean_cov} );
        }
        if ( $pObs <= 0.25 || $pObs >= 0.75 ) { $skip = 1; }
    }

    #if largest cluster is not the closest to 0, then skip as we can't accurately call integer genotypes (at the moment!)
    if ( $max_index != $index || $skip ) {
        for ( my $i = 0 ; $i <= $#sorted ; $i++ ) {
            if    ( $opts{output_format} eq "GL" ) { $$key{ $sorted[$i] }{geno} = "./."; }
            elsif ( $opts{output_format} eq "CN" ) { $$key{ $sorted[$i] }{geno} = "."; }
        }
    }
    else {
        if    ( $opts{output_format} eq "GL" ) { $$key{ $sorted[$index] }{geno} = $genos[$zero]; }
        elsif ( $opts{output_format} eq "CN" ) { $$key{ $sorted[$index] }{geno} = $cn; }
        warn "\tcluster $$key{$sorted[$index]}{mean} assigned to 0/0\n" if $opts{verbose};
        for ( my $i = $index + 1 ; $i <= $#sorted ; $i++ ) {
            $zero++;
            $cn++;
            if ( $opts{output_format} eq "GL" ) {
                if ( defined( $genos[$zero] ) ) {
                    $$key{ $sorted[$i] }{geno} = $genos[$zero];
                }
                else {
                }
                $$key{ $sorted[$i] }{geno} = "./.";
                warn "\tcluster $$key{$sorted[$i]}{mean} assigned to $genos[$zero]\n" if $opts{verbose};
            }
            elsif ( $opts{output_format} eq "CN" ) {
                $$key{ $sorted[$i] }{geno} = $cn;
                warn "\tcluster $$key{$sorted[$i]}{mean} assigned to $cn\n" if $opts{verbose};
            }
        }
        $zero = 3;
        $cn   = 2;
        for ( my $i = $index - 1 ; $i >= 0 ; $i-- ) {
            $zero--;
            $cn--;
            if ( $cn < 0 ) { $cn = 0; }

            if ( $opts{output_format} eq "GL" ) {
                if ( defined( $genos[$zero] ) ) {
                    $$key{ $sorted[$i] }{geno} = $genos[$zero];
                }
                else {
                    $$key{ $sorted[$i] }{geno} = "./.";
                }
                warn "\tcluster $$key{$sorted[$i]}{mean} assigned to $genos[$zero]\n" if $opts{verbose};
            }
            elsif ( $opts{output_format} eq "CN" ) {
                $$key{ $sorted[$i] }{geno} = $cn;
                warn "\tcluster $$key{$sorted[$i]}{mean} assigned to $cn\n" if $opts{verbose};
            }
        }
    }
}

sub round {
    my $float = shift;
    return int( $float + $float / abs( $float * 2 ) );
}

sub factorial {
    my $n = shift;
    my $t = 1;

    $t *= $_ foreach 1 .. $n;
    return $t;
}

sub poisson {
    my $k      = shift;
    my $lambda = shift;
    my $e      = 2.71828;
    return ( ( $lambda**$k ) * $e**( -1 * $lambda ) ) / factorial($k);
}

sub getMeanCov {
    my ( $chr, $pos, $end ) = @_;

    my $num_reads = `$opts{samtools} view -q $opts{min_mapq} $opts{cgh_ref_seq} $chr:$pos-$end | wc -l`;
    my $len       = $end - $pos + 1;
    my $coverage  = $num_reads * $opts{read_len} / $len;
    if   ( $coverage == 0 ) { return 0; }
    else                    { return round($coverage); }
}

sub median {
    my @vals = sort { $a <=> $b } @_;
    my $len = @vals;
    if ( $len % 2 ) {    #odd?
        return $vals[ int( $len / 2 ) ];
    }
    else {               #even
        return ( $vals[ int( $len / 2 ) - 1 ] + $vals[ int( $len / 2 ) ] ) / 2;
    }
}

sub stdev {
    my @vals    = @_;
    my $mean    = mean(@vals);
    my $sqtotal = 0;
    foreach my $val (@vals) {
        $sqtotal += ( $mean - $val )**2;
    }
    my $std = ( $sqtotal / ( @vals - 1 ) )**0.5;
    return $std;
}

sub mean {
    my @vals  = @_;
    my $total = 0;
    foreach my $val (@vals) {
        $total += $val;
    }
    return $total / @vals;
}

sub format_line_hash {

    #modified from Vcf.pm
    my ( $record, $ac, $an, $numProbes ) = @_;

    # CHROM  POS     ID      REF
    my $out;
    $out .= $$record{CHROM} . "\t";
    $out .= $$record{POS} . "\t";
    $out .= ( defined $$record{ID} ? $$record{ID} : '.' ) . "\t";
    $out .= $$record{REF} . "\t";

    # ALT
    $out .= join( ',', @{ $$record{ALT} } ? @{ $$record{ALT} } : '.' );

    # QUAL
    $out .= "\t" . $$record{QUAL};

    # FILTER
    $out .= "\t" . join( ';', $$record{FILTER} ? @{ $$record{FILTER} } : '.' );

    my @info;
    while ( my ( $key, $value ) = each %{ $$record{INFO} } ) {
        if ( defined $value ) {
            if ( $key eq 'AN' ) {
                $value = $an;
            }
            elsif ( $key eq 'AC' ) {
                $value = $ac;
            }
            elsif ( $key eq 'AF' ) {
                if ( $an > 0 ) {
                    $value = sprintf( "%0.4f", $ac / $an );
                }
                else {
                    $value = ".";
                }
            }
            elsif ( $key =~ /_AF/ ) { next; }    #don't have data for population at this point so skip
            push @info, "$key=$value";
        }
        elsif ( $key ne '.' ) {
            push @info, $key;
        }
    }
    push @info, "NPROBES=$numProbes";
    if ( !@info ) { push @info, '.'; }
    $out .= "\t" . join( ';', sort @info );

    # FORMAT
    $out .= "\t$opts{output_format}";
    return $out;
}

sub usage {
    my $version = shift;
    printf("\n");
    printf( "%-9s %s\n", "Program:", "canny.pl" );
    printf( "%-9s %s\n", "Version:", "$version" );
    printf("\n");
    printf( "%-9s %s\n", "Usage:", "canny.pl [options]" );
    printf("\n");
    printf( "%-9s %-35s %s\n", "Options:", "--input_filename=[filename]",         "Input alignment file in VCF format (required if input_range not set)" );
    printf( "%-9s %-35s %s\n", "",         "--input_range=[chr:start-end]",       "Input query range (optional)" );
    printf( "%-9s %-35s %s\n", "",         "--output_filename=[filename]",        "Output file (default stdout)" );
    printf( "%-9s %-35s %s\n", "",         "--ratio_filename=[filename]",         "Probe file contatining log2 ratios per sample, indexed with tabix (required)" );
    printf( "%-9s %-35s %s\n", "",         "--sample_filename=[filename]",        "List of samples in same order as ratio_filename (required)" );
    printf( "%-9s %-35s %s\n", "",         "--work_dir=[directory]",              "Working directory (default .)" );
    printf( "%-9s %-35s %s\n", "",         "--bandwidth=[float]",                 "Bandwidth value for mean shift clustering (default 0.08)" );
    printf( "%-9s %-35s %s\n", "",         "--threshold=[float]",                 "Threshold for merging neighboring clusters (default 0.3)" );
    printf( "%-9s %-35s %s\n", "",         "--cluster_label=[cluster,closest]",   "Label points based on cluster assignment or closest mean (default cluster)" );
    printf( "%-9s %-35s %s\n", "",         "--min_num_probes=[integer]",          "Minimum number of probes to consider site for genotyping (default 5)" );
    printf( "%-9s %-35s %s\n", "",         "--max_cluster_std=[float]",           "Filter sites with largest cluster exceeding standard deviation (default 0.5)" );
    printf( "%-9s %-35s %s\n", "",         "--use_existing_clusters=[directory]", "Use previously determined clusters to save computational time (optional)" );
    printf( "%-9s %-35s %s\n", "",         "--cgh_ref_seq=[filename]",            "Use CGH reference genome to assess sequence depth and correct cluster assignment (optional)" );
    printf( "%-9s %-35s %s\n", "",         "--cgh_ref_mean_cov=[float]",          "Mean coverage for CGH reference genome (required if --cgh_reference_sequence used)" );
    printf( "%-9s %-35s %s\n", "",         "--output_format=[GL,CN],",            "Output format for genotypes (default CN)" );
    printf( "%-9s %-35s %s\n", "",         "--plot_results",                      "Create plots for histograms with fitted densities" );
    printf("\n");
}

sub checkOptions {
    my $optResult = shift;
    my $opts      = shift;
    my $version   = shift;

    if ( !$optResult || $$opts{help} ) {
        usage($version);
        exit;
    }

    if ( !defined( $$opts{input_filename} ) && !defined( $$opts{input_range} ) ) {
        print "\n***ERROR***\t--input_filename or --input_range is required\n";
        usage($version);
        exit;
    }
    elsif ( defined( $$opts{input_filename} ) && defined( $$opts{input_range} ) ) {
        print "\n***ERROR***\t--input_filename and --input_range are mutually exclusive\n";
        usage($version);
        exit;
    }
    elsif ( defined( $$opts{input_range} ) ) {
        my ( $chr, $pos, $end ) = $$opts{input_range} =~ /(\S+):(\d+)-(\d+)/;
        if ( !defined($chr) || !defined($pos) || !defined($end) ) {
            print "\n***ERROR***\t$$opts{input_range} appears malformed\n";
            usage($version);
            exit;
        }
    }
    elsif ( defined( $$opts{cgh_ref_seq} ) && !-e ( $$opts{cgh_ref_seq} ) ) {
        print "\n***ERROR***\t--cgh_ref_seq filename does not exist\n";
        usage($version);
        exit;
    }
    elsif ( defined( $$opts{cgh_ref_seq} ) && !defined( $$opts{cgh_ref_mean_cov} ) ) {
        print "\n***ERROR***\t--cgh_ref_mean_cov is required when --cgh_ref_seq is used\n";
        usage($version);
        exit;
    }
}
