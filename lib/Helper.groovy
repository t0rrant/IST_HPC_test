class Help {

    static def start_info(Map info, String time, String profile) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.4.0"
        println ""
        if (info.containsKey("fastq")){
        int nsamples = info.fastq / 2
        println " Input FastQ                 : $info.fastq"
        println " Input samples               : $nsamples"
        }
        if (info.containsKey("fasta")){
        println " Input Fasta                 : $info.fasta"
        }
        if (info.containsKey("accessions")){
        println " Input accessions            : $info.accessions"
        }
        println " Reports are found in        : ./reports"
        println " Results are found in        : ./results"
        println " Profile                     : $profile"
        println ""
        println "Starting pipeline at $time"
        println ""

    }

    static void complete_info(nextflow.script.WorkflowMetadata wf) {

        println ""
        println "Pipeline execution summary"
        println "=========================="
        println "Completed at                 : $wf.complete"
        println "Duration                     : $wf.duration"
        println "Success                      : $wf.success"
        println "Work directory               : $wf.workDir"
        println "Exit status                  : $wf.exitStatus"
        println ""

    }

    static def print_help(Map params) {

        println ""
        println "============================================================"
        println "                F L O W C R A F T"
        println "============================================================"
        println "Built using flowcraft v1.4.0"
        println ""
        println ""
        println "Usage: "
        println "    nextflow run simple_innuca.nf"
        println ""
        println "       --accessions                Path file with accessions, one perline. (default: $params.fastq) (default: null)"
        println "       "
        println "       Component 'READS_DOWNLOAD_1_1'"
        println "       ------------------------------"
        println "       --asperaKey_1_1             Downloads fastq accessions from ENA using Aspera Connect by providing the private-key file 'asperaweb_id_dsa.openssh' normally found in ~/.aspera/connect/etc/asperaweb_id_dsa.openssh  (default: null)"
        println "       "
        println "       Component 'DOWNSAMPLE_FASTQ_1_2'"
        println "       --------------------------------"
        println "       --genomeSize_1_2            Genome size estimate for the samples in Mb. It is used to estimate the coverage (default: 1)"
        println "       --depth_1_2                 Maximum estimated depth coverage allowed. FastQ with higher estimated depth will be subsampled to this value. (default: 100)"
        println "       --seed_1_2                  The seed number for seqtk. By default it is 100and should be equal for both pairs of reads. (default: 100)"
        println "       --clearInput_1_2            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       "
        println "       Component 'INTEGRITY_COVERAGE_1_3'"
        println "       ----------------------------------"
        println "       --genomeSize_1_3            Genome size estimate for the samples in Mb. It is used to estimate the coverage and other assembly parameters andchecks (default: 1)"
        println "       --minCoverage_1_3           Minimum coverage for a sample to proceed. By default it's setto 0 to allow any coverage (default: 0)"
        println "       "
        println "       Component 'FASTQC_TRIMMOMATIC_1_4'"
        println "       ----------------------------------"
        println "       --adapters_1_4              Path to adapters files, if any. (default: 'None')"
        println "       --trimSlidingWindow_1_4     Perform sliding window trimming, cutting once the average quality within the window falls below a threshold. (default: '5:20')"
        println "       --trimLeading_1_4           Cut bases off the start of a read, if below a threshold quality. (default: 3)"
        println "       --trimTrailing_1_4          Cut bases of the end of a read, if below a threshold quality. (default: 3)"
        println "       --trimMinLength_1_4         Drop the read if it is below a specified length. (default: 55)"
        println "       --clearInput_1_4            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       "
        println "       Component 'CHECK_COVERAGE_1_5'"
        println "       ------------------------------"
        println "       --genomeSize_1_5            Genome size estimate for the samples. It is used to estimate the coverage and other assembly parameters andchecks (default: 2.1)"
        println "       --minCoverage_1_5           Minimum coverage for a sample to proceed. Can be set to0 to allow any coverage (default: 15)"
        println "       "
        println "       Component 'FASTQC_1_6'"
        println "       ----------------------"
        println "       --adapters_1_6              Path to adapters files, if any. (default: 'None')"
        println "       "
        println "       Component 'SPADES_1_7'"
        println "       ----------------------"
        println "       --spadesMinCoverage_1_7     The minimum number of reads to consider an edge in the de Bruijn graph during the assembly (default: 2)"
        println "       --spadesMinKmerCoverage_1_7 Minimum contigs K-mer coverage. After assembly only keep contigs with reported k-mer coverage equal or above this value (default: 2)"
        println "       --spadesKmers_1_7           If 'auto' the SPAdes k-mer lengths will be determined from the maximum read length of each assembly. If 'default', SPAdes will use the default k-mer lengths.  (default: 'auto')"
        println "       --clearInput_1_7            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       --disableRR_1_7             disables repeat resolution stage of assembling. (default: false)"
        println "       "
        println "       Component 'ASSEMBLY_MAPPING_1_8'"
        println "       --------------------------------"
        println "       --minAssemblyCoverage_1_8   In auto, the default minimum coverage for each assembled contig is 1/3 of the assembly mean coverage or 10x, if the mean coverage is below 10x (default: 'auto')"
        println "       --AMaxContigs_1_8           A warning is issued if the number of contigs is overthis threshold. (default: 100)"
        println "       --genomeSize_1_8            Genome size estimate for the samples. It is used to check the ratio of contig number per genome MB (default: 2.1)"
        println "       "
        println "       Component 'PILON_1_9'"
        println "       ---------------------"
        println "       --clearInput_1_9            Permanently removes temporary input files. This option is only useful to remove temporary files in large workflows and prevents nextflow's resume functionality. Use with caution. (default: false)"
        println "       "
        println "       Component 'MLST_1_10'"
        println "       ---------------------"
        println "       --mlstSpecies_1_10          Specify the expected species for MLST checking. (default: null)"
        
    }

}

class CollectInitialMetadata {

    public static void print_metadata(nextflow.script.WorkflowMetadata workflow){

        def treeDag = new File(".treeDag.json").text
        def forkTree = new File(".forkTree.json").text

        def metadataJson = "{'nfMetadata':{'scriptId':'${workflow.scriptId}',\
'scriptName':'${workflow.scriptName}',\
'profile':'${workflow.profile}',\
'container':'${workflow.container}',\
'containerEngine':'${workflow.containerEngine}',\
'commandLine':'${workflow.commandLine}',\
'runName':'${workflow.runName}',\
'sessionId':'${workflow.sessionId}',\
'projectDir':'${workflow.projectDir}',\
'launchDir':'${workflow.launchDir}',\
'startTime':'${workflow.start}',\
'dag':${treeDag},\
'forks':${forkTree}}}"

        def json = metadataJson.replaceAll("'", '"')

        def jsonFile = new File(".metadata.json")
        jsonFile.write json
    }
}