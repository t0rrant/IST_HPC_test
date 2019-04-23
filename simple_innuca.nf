#!/usr/bin/env nextflow

import Helper
import CollectInitialMetadata

// Pipeline version
if (workflow.commitId){
    version = "0.1 $workflow.revision"
} else {
    version = "0.1 (local version)"
}

params.help = false
if (params.help){
    Help.print_help(params)
    exit 0
}

def infoMap = [:]
if (params.containsKey("fastq")){
    infoMap.put("fastq", file(params.fastq).size())
}
if (params.containsKey("fasta")){
    if (file(params.fasta) instanceof LinkedList){
        infoMap.put("fasta", file(params.fasta).size())
    } else {
        infoMap.put("fasta", 1) 
    }
}
if (params.containsKey("accessions")){
    // checks if params.accessions is different from null
    if (params.accessions) {
        BufferedReader reader = new BufferedReader(new FileReader(params.accessions));
        int lines = 0;
        while (reader.readLine() != null) lines++;
        reader.close();
        infoMap.put("accessions", lines)
    }
}

Help.start_info(infoMap, "$workflow.start", "$workflow.profile")
CollectInitialMetadata.print_metadata(workflow)
    

// Placeholder for main input channels
if (!params.accessions){ exit 1, "'accessions' parameter missing" }

IN_accessions_raw = Channel.fromPath(params.accessions).ifEmpty { exit 1, "No accessions file provided with path:'${params.accessions}'" }

// Placeholder for secondary input channels


// Placeholder for extra input channels


// Placeholder to fork the raw input channel

IN_accessions_raw.set{ reads_download_in_1_0 }


if (params.asperaKey_1_1){
    if (file(params.asperaKey_1_1).exists()){
        IN_asperaKey_1_1 = Channel.fromPath(params.asperaKey_1_1)
    } else {
        IN_asperaKey_1_1 = Channel.value("")
    }
} else {
    IN_asperaKey_1_1 = Channel.value("")
}

process reads_download_1_1 {

        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_1 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_1 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_1 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId reads_download_1_1 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { accession_id }
    publishDir "reads", pattern: "${accession_id}/*fq.gz"
    maxRetries 1

    input:
    set val(accession_id), val(name) from reads_download_in_1_0.splitText(){ it.trim() }.unique().filter{ it != "" }.map{ it.split().length > 1 ? ["accession": it.split()[0], "name": it.split()[1]] : [it.split()[0], null] }
    each file(aspera_key) from IN_asperaKey_1_1

    output:
    set val({ "$name" != "null" ? "$name" : "$accession_id" }), file("${accession_id}/*fq.gz") optional true into reads_download_out_1_0
    set accession_id, val("1_1_reads_download"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_reads_download_1_1
set accession_id, val("reads_download_1_1"), val("1_1"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_reads_download_1_1
file ".versions"

    script:
    """
    {
        # getSeqENA requires accession numbers to be provided as a text file
        echo "${accession_id}" >> accession_file.txt
        # Set default status value. It will be overwritten if anything goes wrong
        echo "pass" > ".status"

        if [ -f $aspera_key ]; then
            asperaOpt="-a $aspera_key"
        else
            asperaOpt=""
        fi

        getSeqENA.py -l accession_file.txt \$asperaOpt -o ./ --SRAopt --downloadCramBam

        # If a name has been provided along with the accession, rename the
        # fastq files.
        if [ $name != null ];
        then
            echo renaming pattern '${accession_id}' to '${name}' && cd ${accession_id} && rename "s/${accession_id}/${name}/" *.gz
        fi
    } || {
        # If exit code other than 0
        if [ \$? -eq 0 ]
        then
            echo "pass" > .status
        else
            echo "fail" > .status
            echo "Could not download accession $accession_id" > .fail
        fi
    }
    version_str="{'version':[{'program':'getSeqENA.py','version':'1.3'}]}"
    echo \$version_str > .versions
    """

}



IN_genome_size_1_2 = Channel.value(params.genomeSize_1_2)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize_1_2}'")}

IN_depth_1_2 = Channel.value(params.depth_1_2)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The depth parameter must be a number or a float. Provided value: '${params.depth_1_2}'")}

IN_seed_1_2 = Channel.value(params.seed_1_2)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The seed parameter must be a number or a float. Provided value: '${params.seed_1_2}'")}

clear = params.clearInput_1_2 ? "true" : "false"
checkpointClear_1_2 = Channel.value(clear)

process downsample_fastq_1_2 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_2 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_2 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId downsample_fastq_1_2 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { "${sample_id}" }
    publishDir "results/downsample_fastq_1_2/", pattern: "_ss.*"

    input:
    set sample_id, file(fastq_pair) from reads_download_out_1_0
    val gsize from IN_genome_size_1_2
    val depth from IN_depth_1_2
    val seed from IN_seed_1_2
    val clear from checkpointClear_1_2

    output:
    set sample_id, file('*_ss.*') into downsample_fastq_out_1_1
    set sample_id, val("1_2_downsample_fastq"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_downsample_fastq_1_2
set sample_id, val("downsample_fastq_1_2"), val("1_2"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_downsample_fastq_1_2
file ".versions"

    script:
    template "downsample_fastq.py"

}


IN_genome_size_1_3 = Channel.value(params.genomeSize_1_3)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize__1_3}'")}

IN_min_coverage_1_3 = Channel.value(params.minCoverage_1_3)
    .map{it -> it.toString().isNumber() ? it : exit(1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage__1_3}'")}

process integrity_coverage_1_3 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_3 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_3 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId integrity_coverage_1_3 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(fastq_pair) from downsample_fastq_out_1_1
    val gsize from IN_genome_size_1_3
    val cov from IN_min_coverage_1_3
    // This channel is for the custom options of the integrity_coverage.py
    // script. See the script's documentation for more information.
    val opts from Channel.value('')

    output:
    set sample_id,
        file(fastq_pair),
        file('*_encoding'),
        file('*_phred'),
        file('*_coverage'),
        file('*_max_len') into MAIN_integrity_1_3
    file('*_report') optional true into LOG_report_coverage1_1_3
    set sample_id, val("1_3_integrity_coverage"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_integrity_coverage_1_3
set sample_id, val("integrity_coverage_1_3"), val("1_3"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_integrity_coverage_1_3
file ".versions"

    script:
    template "integrity_coverage.py"

}

// TRIAGE OF CORRUPTED SAMPLES
LOG_corrupted_1_3 = Channel.create()
MAIN_PreCoverageCheck_1_3 = Channel.create()
// Corrupted samples have the 2nd value with 'corrupt'
MAIN_integrity_1_3.choice(LOG_corrupted_1_3, MAIN_PreCoverageCheck_1_3) {
    a -> a[2].text == "corrupt" ? 0 : 1
}

// TRIAGE OF LOW COVERAGE SAMPLES
integrity_coverage_out_1_2 = Channel.create()
SIDE_phred_1_3 = Channel.create()
SIDE_max_len_1_3 = Channel.create()

MAIN_PreCoverageCheck_1_3
// Low coverage samples have the 4th value of the Channel with 'fail'
    .filter{ it[4].text != "fail" }
// For the channel to proceed with FastQ in 'sample_good' and the
// Phred scores for each sample in 'SIDE_phred'
    .separate(integrity_coverage_out_1_2, SIDE_phred_1_3, SIDE_max_len_1_3){
        a -> [ [a[0], a[1]], [a[0], a[3].text], [a[0], a[5].text]  ]
    }

/** REPORT_COVERAGE - PLUG-IN
This process will report the expected coverage for each non-corrupted sample
and write the results to 'reports/coverage/estimated_coverage_initial.csv'
*/
process report_coverage_1_3 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage_1_3/'

    input:
    file(report) from LOG_report_coverage1_1_3.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_initial.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_initial.csv
    cat $report >> estimated_coverage_initial.csv
    """
}

/** REPORT_CORRUPT - PLUG-IN
This process will report the corrupted samples and write the results to
'reports/corrupted/corrupted_samples.txt'
*/
process report_corrupt_1_3 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/corrupted_1_3/'

    input:
    val sample_id from LOG_corrupted_1_3.collect{it[0]}

    output:
    file 'corrupted_samples.txt'

    """
    echo ${sample_id.join(",")} | tr "," "\n" >> corrupted_samples.txt
    """

}


SIDE_phred_1_3.set{ SIDE_phred_1_4 }


// Check sliding window parameter
if ( params.trimSlidingWindow_1_4.toString().split(":").size() != 2 ){
    exit 1, "'trimSlidingWindow_1_4' parameter must contain two values separated by a ':'. Provided value: '${params.trimSlidingWindow_1_4}'"
}
if ( !params.trimLeading_1_4.toString().isNumber() ){
    exit 1, "'trimLeading_1_4' parameter must be a number. Provide value: '${params.trimLeading_1_4}'"
}
if ( !params.trimTrailing_1_4.toString().isNumber() ){
    exit 1, "'trimTrailing_1_4' parameter must be a number. Provide value: '${params.trimTrailing_1_4}'"
}
if ( !params.trimMinLength_1_4.toString().isNumber() ){
    exit 1, "'trimMinLength_1_4' parameter must be a number. Provide value: '${params.trimMinLength_1_4}'"
}

IN_trimmomatic_opts_1_4 = Channel.value([params.trimSlidingWindow_1_4,params.trimLeading_1_4,params.trimTrailing_1_4,params.trimMinLength_1_4])
IN_adapters_1_4 = Channel.value(params.adapters_1_4)

clear = params.clearInput_1_4 ? "true" : "false"
checkpointClear_1_4 = Channel.value(clear)

process fastqc_1_4 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_4 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir "reports/fastqc_1_4/", pattern: "*.html"

    input:
    set sample_id, file(fastq_pair) from integrity_coverage_out_1_2
    val ad from Channel.value('None')

    output:
    set sample_id, file(fastq_pair), file('pair_1*'), file('pair_2*') into MAIN_fastqc_out_1_4
    file "*html"
    set sample_id, val("1_4_fastqc"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc_1_4
set sample_id, val("fastqc_1_4"), val("1_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc_1_4
file ".versions"

    script:
    template "fastqc.py"
}

/** FASTQC_REPORT - MAIN
This process will parse the result files from a FastQC analyses and output
the optimal_trim information for Trimmomatic
*/
process fastqc_report_1_4 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_4 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc_1_4/run_1/', pattern: '*summary.txt', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_1_4
    val opts from Channel.value("--ignore-tests")

    output:
    set sample_id, file(fastq_pair), 'optimal_trim', ".status" into _MAIN_fastqc_trim_1_4
    file '*_trim_report' into LOG_trim_1_4
    file "*_status_report" into LOG_fastqc_report_1_4
    file "${sample_id}_*_summary.txt" optional true
    set sample_id, val("1_4_fastqc_report"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc_report_1_4
set sample_id, val("fastqc_report_1_4"), val("1_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc_report_1_4
file ".versions"

    script:
    template "fastqc_report.py"

}

MAIN_fastqc_trim_1_4 = Channel.create()
_MAIN_fastqc_trim_1_4
        .filter{ it[3].text == "pass" }
        .map{ [it[0], it[1], file(it[2]).text] }
        .into(MAIN_fastqc_trim_1_4)


/** TRIM_REPORT - PLUG-IN
This will collect the optimal trim points assessed by the fastqc_report
process and write the results of all samples in a single csv file
*/
process trim_report_1_4 {

    publishDir 'reports/fastqc_1_4/', mode: 'copy'

    input:
    file trim from LOG_trim_1_4.collect()

    output:
    file "FastQC_trim_report.csv"

    """
    echo Sample,Trim begin, Trim end >> FastQC_trim_report.csv
    cat $trim >> FastQC_trim_report.csv
    """
}


process compile_fastqc_status_1_4 {

    publishDir 'reports/fastqc_1_4/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_1_4.collect()

    output:
    file 'FastQC_1run_report.csv'

    """
    echo Sample, Failed? >> FastQC_1run_report.csv
    cat $rep >> FastQC_1run_report.csv
    """

}


/** TRIMMOMATIC - MAIN
This process will execute trimmomatic. Currently, the main channel requires
information on the trim_range and phred score.
*/
process trimmomatic_1_4 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_4 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_4 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_trimmomatic_1_4 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    publishDir "results/trimmomatic_1_4", pattern: "*.gz"

    input:
    set sample_id, file(fastq_pair), trim_range, phred from MAIN_fastqc_trim_1_4.join(SIDE_phred_1_4)
    val opts from IN_trimmomatic_opts_1_4
    val ad from IN_adapters_1_4
    val clear from checkpointClear_1_4

    output:
    set sample_id, "${sample_id}_*trim.fastq.gz" into fastqc_trimmomatic_out_1_3
    file 'trimmomatic_report.csv'
    set sample_id, val("1_4_trimmomatic"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_trimmomatic_1_4
set sample_id, val("trimmomatic_1_4"), val("1_4"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_trimmomatic_1_4
file ".versions"

    script:
    template "trimmomatic.py"

}



IN_genome_size_1_5 = Channel.value(params.genomeSize_1_5)
    .map{it -> it.toString().isNumber() ? it : exit (1, "The genomeSize parameter must be a number or a float. Provided value: '${params.genomeSize_1_5}'")}
IN_min_coverage_1_5 = Channel.value(params.minCoverage_1_5)
    .map{it -> it.toString().isNumber() ? it : exit (1, "The minCoverage parameter must be a number or a float. Provided value: '${params.minCoverage_1_5}'")}

process integrity_coverage2_1_5 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_5 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_5 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId check_coverage_1_5 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    cpus 1

    input:
    set sample_id, file(fastq_pair) from fastqc_trimmomatic_out_1_3
    val gsize from IN_genome_size_1_5
    val cov from IN_min_coverage_1_5
    // Use -e option for skipping encoding guess
    val opts from Channel.value('-e')

    output:
    set sample_id,
        file(fastq_pair),
        file('*_coverage'),
        file('*_max_len') optional true into MAIN_integrity_1_5
    file('*_report') into LOG_report_coverage_1_5
    set sample_id, val("1_5_check_coverage"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_check_coverage_1_5
set sample_id, val("check_coverage_1_5"), val("1_5"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_check_coverage_1_5
file ".versions"

    script:
    template "integrity_coverage.py"
}

check_coverage_out_1_4 = Channel.create()
SIDE_max_len_1_5 = Channel.create()

MAIN_integrity_1_5
    .filter{ it[2].text != "fail" }
    .separate(check_coverage_out_1_4, SIDE_max_len_1_5){
        a -> [ [a[0], a[1]], [a[0], a[3].text]]
    }


process report_coverage2_1_5 {

    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/coverage_1_5/'

    input:
    file(report) from LOG_report_coverage_1_5.filter{ it.text != "corrupt" }.collect()

    output:
    file 'estimated_coverage_second.csv'

    """
    echo Sample,Estimated coverage,Status >> estimated_coverage_second.csv
    cat $report >> estimated_coverage_second.csv
    """
}


SIDE_max_len_1_5.set{ SIDE_max_len_1_7 }


IN_adapters_1_6 = Channel.value(params.adapters_1_6)

process fastqc2_1_6 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_6 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_6 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_6 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_1_6 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir "reports/fastqc_1_6/", pattern: "*.html"

    input:
    set sample_id, file(fastq_pair) from check_coverage_out_1_4
    val ad from IN_adapters_1_6

    output:
    set sample_id, file(fastq_pair), file('pair_1*'), file('pair_2*') into MAIN_fastqc_out_1_6
    file "*html"
    set sample_id, val("1_6_fastqc2"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc2_1_6
set sample_id, val("fastqc2_1_6"), val("1_6"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc2_1_6
file ".versions"

    script:
    template "fastqc.py"
}


process fastqc2_report_1_6 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_6 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_6 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_6 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId fastqc_1_6 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1
    publishDir 'reports/fastqc_1_6/run_2/', pattern: '*summary.txt', mode: 'copy'

    input:
    set sample_id, file(fastq_pair), file(result_p1), file(result_p2) from MAIN_fastqc_out_1_6
    val opts from Channel.value("")

    output:
    set sample_id, file(fastq_pair), '.status' into MAIN_fastqc_report_1_6
    file "*_status_report" into LOG_fastqc_report_1_6
    file "${sample_id}_*_summary.txt" optional true
    set sample_id, val("1_6_fastqc2_report"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_fastqc2_report_1_6
set sample_id, val("fastqc2_report_1_6"), val("1_6"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_fastqc2_report_1_6
file ".versions"

    script:
    template "fastqc_report.py"

}


process compile_fastqc_status2_1_6 {

    publishDir 'reports/fastqc_1_6/', mode: 'copy'

    input:
    file rep from LOG_fastqc_report_1_6.collect()

    output:
    file 'FastQC_2run_report.csv'

    """
    echo Sample, Failed? >> FastQC_2run_report.csv
    cat $rep >> FastQC_2run_report.csv
    """

}

_fastqc_out_1_5 = Channel.create()

MAIN_fastqc_report_1_6
        .filter{ it[2].text == "pass" }
        .map{ [it[0], it[1]] }
        .into(_fastqc_out_1_5)


_fastqc_out_1_5.into{ fastqc_out_1_5;_LAST_fastq_1_8 }


if ( !params.spadesMinCoverage_1_7.toString().isNumber() ){
    exit 1, "'spadesMinCoverage_1_7' parameter must be a number. Provided value: '${params.spadesMinCoverage_1_7}'"
}
if ( !params.spadesMinKmerCoverage_1_7.toString().isNumber()){
    exit 1, "'spadesMinKmerCoverage_1_7' parameter must be a number. Provided value: '${params.spadesMinKmerCoverage_1_7}'"
}

IN_spades_opts_1_7 = Channel.value(
    [params.spadesMinCoverage_1_7,
     params.spadesMinKmerCoverage_1_7
     ])

if ( params.spadesKmers_1_7.toString().split(" ").size() <= 1 ){
    if (params.spadesKmers_1_7.toString() != 'auto'){
        exit 1, "'spadesKmers_1_7' parameter must be a sequence of space separated numbers or 'auto'. Provided value: ${params.spadesKmers_1_7}"
    }
}
IN_spades_kmers_1_7 = Channel.value(params.spadesKmers_1_7)

clear = params.clearInput_1_7 ? "true" : "false"
disable_rr = params.disableRR_1_7 ? "true" : "false"

checkpointClear_1_7 = Channel.value(clear)
disableRR_1_7 = Channel.value(disable_rr)

process spades_1_7 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_7 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_7 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_7 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId spades_1_7 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    publishDir 'results/assembly/spades_1_7/', pattern: '*_spades*.fasta', mode: 'copy'
    publishDir "results/assembly/spades_1_7/$sample_id", pattern: "*.gfa", mode: "copy"
    publishDir "results/assembly/spades_1_7/$sample_id", pattern: "*.fastg", mode: "copy"

    input:
    set sample_id, file(fastq_pair), max_len from fastqc_out_1_5.join(SIDE_max_len_1_7)
    val opts from IN_spades_opts_1_7
    val kmers from IN_spades_kmers_1_7
    val clear from checkpointClear_1_7
    val disable_rr from disableRR_1_7

    output:
    set sample_id, file('*_spades*.fasta') into spades_out_1_6
    file "*.fastg" optional true
    file "*.gfa" into gfa1_1_7
    set sample_id, val("1_7_spades"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_spades_1_7
set sample_id, val("spades_1_7"), val("1_7"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_spades_1_7
file ".versions"

    script:
    template "spades.py"

}



if ( !params.minAssemblyCoverage_1_8.toString().isNumber() ){
    if (params.minAssemblyCoverage_1_8.toString() != 'auto'){
        exit 1, "'minAssemblyCoverage_1_8' parameter must be a number or 'auto'. Provided value: ${params.minAssemblyCoverage_1_8}"
    }
}
if ( !params.AMaxContigs_1_8.toString().isNumber() ){
    exit 1, "'AMaxContigs_1_8' parameter must be a number. Provide value: '${params.AMaxContigs_1_8}'"
}

IN_assembly_mapping_opts_1_8 = Channel.value([params.minAssemblyCoverage_1_8,params.AMaxContigs_1_8])
IN_genome_size_1_8 = Channel.value(params.genomeSize_1_8)


process assembly_mapping_1_8 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_8 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_8 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_8 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId assembly_mapping_1_8 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }

    input:
    set sample_id, file(assembly), file(fastq) from spades_out_1_6.join(_LAST_fastq_1_8)

    output:
    set sample_id, file(assembly), 'coverages.tsv', 'coverage_per_bp.tsv', 'sorted.bam', 'sorted.bam.bai' into MAIN_am_out_1_8
    set sample_id, file("coverage_per_bp.tsv") optional true into SIDE_BpCoverage_1_8
    set sample_id, val("1_8_assembly_mapping"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_assembly_mapping_1_8
set sample_id, val("assembly_mapping_1_8"), val("1_8"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_assembly_mapping_1_8
file ".versions"

    script:
    """
    {
        echo [DEBUG] BUILDING BOWTIE INDEX FOR ASSEMBLY: $assembly >> .command.log 2>&1
        bowtie2-build --threads ${task.cpus} $assembly genome_index >> .command.log 2>&1
        echo [DEBUG] MAPPING READS FROM $fastq >> .command.log 2>&1
        bowtie2 -q --very-sensitive-local --threads ${task.cpus} -x genome_index -1 ${fastq[0]} -2 ${fastq[1]} -S mapping.sam >> .command.log 2>&1
        echo [DEBUG] CONVERTING AND SORTING SAM TO BAM >> .command.log 2>&1
        samtools sort -o sorted.bam -O bam -@ ${task.cpus} mapping.sam && rm *.sam  >> .command.log 2>&1
        echo [DEBUG] CREATING BAM INDEX >> .command.log 2>&1
        samtools index sorted.bam >> .command.log 2>&1
        echo [DEBUG] ESTIMATING READ DEPTH >> .command.log 2>&1
        parallel -j ${task.cpus} samtools depth -ar {} sorted.bam \\> {}.tab  ::: \$(grep ">" $assembly | cut -c 2- | tr " " "_")
        # Insert 0 coverage count in empty files. See Issue #2
        echo [DEBUG] REMOVING EMPTY FILES  >> .command.log 2>&1
        find . -size 0 -print0 | xargs -0 -I{} sh -c 'echo -e 0"\t"0"\t"0 > "{}"'
        echo [DEBUG] COMPILING COVERAGE REPORT  >> .command.log 2>&1
        parallel -j ${task.cpus} echo -n {.} '"\t"' '&&' cut -f3 {} '|' paste -sd+ '|' bc >> coverages.tsv  ::: *.tab
        cat *.tab > coverage_per_bp.tsv
        rm *.tab
        if [ -f "coverages.tsv" ]
        then
            echo pass > .status
        else
            echo fail > .status
        fi
        echo -n "" > .report.json
        echo -n "" > .versions
    } || {
        echo fail > .status
    }
    """
}


/** PROCESS_ASSEMBLY_MAPPING -  MAIN
Processes the results from the assembly_mapping process and filters the
assembly contigs based on coverage and length thresholds.
*/
process process_assembly_mapping_1_8 {

    // Send POST request to platform
    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_8 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_8 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_8 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId assembly_mapping_1_8 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(assembly), file(coverage), file(coverage_bp), file(bam_file), file(bam_index) from MAIN_am_out_1_8
    val opts from IN_assembly_mapping_opts_1_8
    val gsize from IN_genome_size_1_8

    output:
    set sample_id, '*_filt.fasta', 'filtered.bam', 'filtered.bam.bai' into assembly_mapping_out_1_7
    set sample_id, val("1_8_process_am"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_process_am_1_8
set sample_id, val("process_am_1_8"), val("1_8"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_process_am_1_8
file ".versions"

    script:
    template "process_assembly_mapping.py"

}


SIDE_BpCoverage_1_8.set{ SIDE_BpCoverage_1_9 }



clear = params.clearInput_1_9 ? "true" : "false"
checkpointClear_1_9 = Channel.value(clear)

process pilon_1_9 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_9 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_9 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_9 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId pilon_1_9 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    echo false
    publishDir 'results/assembly/pilon_1_9/', mode: 'copy', pattern: "*.fasta"

    input:
    set sample_id, file(assembly), file(bam_file), file(bam_index) from assembly_mapping_out_1_7
    val clear from checkpointClear_1_9

    output:
    set sample_id, '*_polished.fasta' into pilon_out_1_8, pilon_report_1_9
    set sample_id, val("1_9_pilon"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_pilon_1_9
set sample_id, val("pilon_1_9"), val("1_9"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_pilon_1_9
file ".versions"

    script:
    """
    {
        pilon_mem=${String.valueOf(task.memory).substring(0, String.valueOf(task.memory).length() - 1).replaceAll("\\s", "")}
        java -jar -Xms256m -Xmx\${pilon_mem} /NGStools/pilon-1.22.jar --genome $assembly --frags $bam_file --output ${assembly.name.replaceFirst(~/\.[^\.]+$/, '')}_polished --changes --threads $task.cpus >> .command.log 2>&1
        echo pass > .status

        if [ "$clear" = "true" ];
        then
            work_regex=".*/work/.{2}/.{30}/.*"
            assembly_file=\$(readlink -f \$(pwd)/${assembly})
            bam_file=\$(readlink -f \$(pwd)/${bam_file})
            if [[ "\$assembly_file" =~ \$work_regex ]]; then
                rm \$assembly_file \$bam_file
            fi
        fi

    } || {
        echo fail > .status
    }
    """

}

process pilon_report_1_9 {

    
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh"
        afterScript "report_POST.sh $params.projectId $params.pipelineId 1_9 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId pilon_1_9 \"$params.platformSpecies\" false"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh"
        }
    

    tag { sample_id }

    input:
    set sample_id, file(assembly), file(coverage_bp) from pilon_report_1_9.join(SIDE_BpCoverage_1_9)

    output:
    file "*_assembly_report.csv" into pilon_report_out_1_9
    set sample_id, val("1_9_pilon_report"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_pilon_report_1_9
set sample_id, val("pilon_report_1_9"), val("1_9"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_pilon_report_1_9
file ".versions"

    script:
    template "assembly_report.py"

}


process compile_pilon_report_1_9 {

    publishDir "reports/assembly/pilon_1_9/", mode: 'copy'

    input:
    file(report) from pilon_report_out_1_9.collect()

    output:
    file "pilon_assembly_report.csv"

    """
    echo Sample,Number of contigs,Average contig size,N50,Total assembly length,GC content,Missing data > pilon_assembly_report.csv
    cat $report >> pilon_assembly_report.csv
    """
}




process mlst_1_10 {

    // Send POST request to platform
        if ( params.platformHTTP != null ) {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH; set_dotfiles.sh; startup_POST.sh $params.projectId $params.pipelineId 1_10 $params.platformHTTP"
        afterScript "final_POST.sh $params.projectId $params.pipelineId 1_10 $params.platformHTTP; report_POST.sh $params.projectId $params.pipelineId 1_10 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId mlst_1_10 \"$params.platformSpecies\" true"
    } else {
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; set_dotfiles.sh"
        }

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(assembly) from pilon_out_1_8

    output:
    file '*.mlst.txt' into LOG_mlst_1_10
    set sample_id, file(assembly), file(".status") into MAIN_mlst_out_1_10
    set sample_id, val("1_10_mlst"), file(".status"), file(".warning"), file(".fail"), file(".command.log") into STATUS_mlst_1_10
set sample_id, val("mlst_1_10"), val("1_10"), file(".report.json"), file(".versions"), file(".command.trace") into REPORT_mlst_1_10
file ".versions"

    script:
    """
    {
        expectedSpecies=${params.mlstSpecies_1_10}
        mlst $assembly >> ${sample_id}.mlst.txt
        mlstSpecies=\$(cat *.mlst.txt | cut -f2)
        json_str="{'expectedSpecies':\'\$expectedSpecies\',\
            'species':'\$mlstSpecies',\
            'st':'\$(cat *.mlst.txt | cut -f3)',\
            'tableRow':[{'sample':'${sample_id}','data':[\
                {'header':'MLST species','value':'\$mlstSpecies','table':'typing'},\
                {'header':'MLST ST','value':'\$(cat *.mlst.txt | cut -f3)','table':'typing'}]}]}"
        echo \$json_str > .report.json

        if [ ! \$mlstSpecies = \$expectedSpecies ];
        then
            printf fail > .status
        else
            printf pass > .status
        fi

    } || {
        printf fail > .status
    }
    """
}

process compile_mlst_1_10 {

    publishDir "results/annotation/mlst_1_10/"

    input:
    file res from LOG_mlst_1_10.collect()

    output:
    file "mlst_report.tsv"

    script:
    """
    cat $res >> mlst_report.tsv
    """
}

mlst_out_1_9 = Channel.create()
MAIN_mlst_out_1_10
    .filter{ it[2].text != "fail" }
    .map{ [it[0], it[1]] }
    .set{ mlst_out_1_9 }





/** STATUS
Reports the status of a sample in any given process.
*/
process status {

    tag { sample_id }
    publishDir "pipeline_status/$task_name"

    input:
    set sample_id, task_name, status, warning, fail, file(log) from STATUS_reads_download_1_1.mix(STATUS_downsample_fastq_1_2,STATUS_integrity_coverage_1_3,STATUS_fastqc_1_4,STATUS_fastqc_report_1_4,STATUS_trimmomatic_1_4,STATUS_check_coverage_1_5,STATUS_fastqc2_1_6,STATUS_fastqc2_report_1_6,STATUS_spades_1_7,STATUS_assembly_mapping_1_8,STATUS_process_am_1_8,STATUS_pilon_1_9,STATUS_pilon_report_1_9,STATUS_mlst_1_10)

    output:
    file '*.status' into master_status
    file '*.warning' into master_warning
    file '*.fail' into master_fail
    file '*.log'

    """
    echo $sample_id, $task_name, \$(cat $status) > ${sample_id}_${task_name}.status
    echo $sample_id, $task_name, \$(cat $warning) > ${sample_id}_${task_name}.warning
    echo $sample_id, $task_name, \$(cat $fail) > ${sample_id}_${task_name}.fail
    echo "\$(cat .command.log)" > ${sample_id}_${task_name}.log
    """
}

process compile_status_buffer {

    input:
    file status from master_status.buffer( size: 5000, remainder: true)
    file warning from master_warning.buffer( size: 5000, remainder: true)
    file fail from master_fail.buffer( size: 5000, remainder: true)

    output:
    file 'master_status_*.csv' into compile_status_buffer
    file 'master_warning_*.csv' into compile_warning_buffer
    file 'master_fail_*.csv' into compile_fail_buffer

    """
    cat $status >> master_status_${task.index}.csv
    cat $warning >> master_warning_${task.index}.csv
    cat $fail >> master_fail_${task.index}.csv
    """
}

process compile_status {

    publishDir 'reports/status'

    input:
    file status from compile_status_buffer.collect()
    file warning from compile_warning_buffer.collect()
    file fail from compile_fail_buffer.collect()

    output:
    file "*.csv"

    """
    cat $status >> master_status.csv
    cat $warning >> master_warning.csv
    cat $fail >> master_fail.csv
    """

}


/** Reports
Compiles the reports from every process
*/
process report {

    tag { sample_id }

    input:
    set sample_id,
            task_name,
            pid,
            report_json,
            version_json,
            trace from REPORT_reads_download_1_1.mix(REPORT_downsample_fastq_1_2,REPORT_integrity_coverage_1_3,REPORT_fastqc_1_4,REPORT_fastqc_report_1_4,REPORT_trimmomatic_1_4,REPORT_check_coverage_1_5,REPORT_fastqc2_1_6,REPORT_fastqc2_report_1_6,REPORT_spades_1_7,REPORT_assembly_mapping_1_8,REPORT_process_am_1_8,REPORT_pilon_1_9,REPORT_pilon_report_1_9,REPORT_mlst_1_10)

    output:
    file "*" optional true into master_report

    """
    prepare_reports.py $report_json $version_json $trace $sample_id $task_name 1 $pid $workflow.scriptId $workflow.runName
    """

}


process compile_reports {

    publishDir "pipeline_report/", mode: "copy"

    if ( params.reportHTTP != null ){
        beforeScript "PATH=${workflow.projectDir}/bin:\$PATH; export PATH;"
        afterScript "metadata_POST.sh $params.projectId $params.pipelineId 0 $params.sampleName $params.reportHTTP $params.currentUserName $params.currentUserId 0 \"$params.platformSpecies\""
    }

    input:
    file report from master_report.collect()
    file forks from Channel.fromPath(".forkTree.json")
    file dag from Channel.fromPath(".treeDag.json")
    file js from Channel.fromPath("${workflow.projectDir}/resources/main.js.zip")

    output:
    file "pipeline_report.json"
    file "pipeline_report.html"
    file "src/main.js"

    script:
    template "compile_reports.py"
}

workflow.onComplete {
  // Display complete message
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
}

workflow.onError {
  // Display error message
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
