//Genome specific
params.genome = ''
params.genomes = []
params.bt2_index = params.genome ? params.genomes[ params.genome ].bt2Index ?: false : false
params.blacklist = params.genome ? params.genomes[ params.genome ].blacklist ?: false : false
params.genesList = params.genome ? params.genomes[ params.genome ].genesList ?: false : false



version = 0.5


def helpMessage() {
	log.info """
		=================================================
            C U T & R U N  P I P E L I N E v${version}
        =================================================
		Author: Michael J. Bale (mib4004@med.cornell.edu)
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run michaelbale/JLab-Flow --input 'project/*_R{1,2}.fastq.gz' --genome mouse -profile singularity 
    Mandatory arguments:
      --input                       Path to input data (must be surrounded with quotes)
      --genome                      Name of Genomes reference (current supported: mouse -- mm10, human -- hg38)
	  -profile                      Name of package manager system (available: docker, singularity, conda); for WCM default -- singularity is recommended, but conda works while docker does not. For minimal - use conda.
      
    Options:
      --executorConfig              Path to config file that contains specifics for execution. Will default to WCM SCU-specific parameters. N.B. for local running use --executorConfig conf/minimal.config
	  --singleSample                Specifies that the input is a single sample and will not generate a PCA graph
      --PCATitle                    Title to be included in PCA graph; must be surrounded with quotes
	  --catLanes                    Tells CnRFlow to take input files and concatenate lanes into single fastq files
	  --name                        Project Name; cannot have whitespace characters
	  --addBEDFilesProfile          Path to csv file with info on additional BED files for generating rStart-rEnd profile plots; csv format: rName,BEDPath
	  --addBEDFilesRefPoint         Path to csv file with info on additional BED files for generating rStart +/- region profile plots; csv format: pName,BEDPath,PlusMinus,pLabel
	  --workDir                     Name of folder to output concatenated fastq files to (not used unless --catLanes)
	  --outdir                      Name of folder to output all results to
	  --genomeAssets                Home directory of where genome-specific files are. Defaults to /athena/josefowiczlab/scratch/szj2001/JLab-Flow_Genomes
	  --singleEnd                   Denotes read input are single read instead of paired end sequencing
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}




log.info """\
		=================================================
            C U T & R U N  P I P E L I N E v${version}
        =================================================
		Author: Michael J. Bale (mib4004@med.cornell.edu)
		
		Project ID: ${params.name}
        Genome: ${params.genome}
        Reads: ${params.input}
        Publish Directory: ${params.outdir}
        """
         .stripIndent()

/*
		 
    getSampleID = {
	    (it =~ /(.+)_S\d+_L\d{3}/)[0][1]
	}

	getGroupID = {
	  (it =~ /^(.+)-.+_S\d+_L\d{3}$/))[0][1]
	}
	
	getTargetID = {
	  (it =~ /^.+-(.+)_S\d+_L\d{3}$/))[0][1]
	}
*/	
if(params.catLanes) {

    getSampleID = {
        (it =~ /(.+)_S\d+_L\d{3}/)[0][1]
    }
    
	if(params.singleEnd) {
	
	  print("y u do dis 2 me andrew")
	    Channel
		  .fromPath(params.input)
		  .map { path -> tuple(getSampleID(path.getName()), path)}
		  .groupTuple()
		  .set { inFq_ch }
		
		process catLanesSE {
		  tag "Concatenating lanes into $params.workDir"
		  publishDir "$params.workDir/$sampleID", mode: 'copy', pattern: '*.gz'
		  label 'small_mem'
		  
		  input:
		  tuple val(sampleID), path(r1) from inFq_ch
		  
		  output:
		  tuple val(sampleID), path("${sampleID}_R1_init.fq.gz") into reads_ch
		  
		  script:
		  """
		  zcat $r1 > ${sampleID}_R1_init.fq
		  gzip  ${sampleID}_R1_init.fq
		  """
		}
	} else {
	  Channel
	  .fromFilePairs(params.input, flat: true)
	  .map { prefix, r1, r2 -> tuple(getSampleID(prefix), r1, r2) }
	  .groupTuple()
	  .set {inFq_ch}
	

      process catLanes {
	    tag "Concatenating lanes into $params.workDir"
		publishDir "$params.workDir/$sampleID", mode: 'copy', pattern: "*.gz"
		label 'small_mem'
		
		input:
		tuple val(sampleID), path(R1), path(R2) from inFq_ch
				
		output:
		tuple val(sampleID), path("${sampleID}_*_init.fq.gz") into reads_ch
		
		script:
		"""
		zcat $R1 > ${sampleID}_R1_init.fq
		zcat $R2 > ${sampleID}_R2_init.fq
		
		gzip ${sampleID}_R1_init.fq
		gzip ${sampleID}_R2_init.fq
		"""
	  }
	}
} else {
	
	if(params.singleEnd) {
	    getSampleID = {
	    (it =~ /(.+)_.+.gz/)[0][1]
		}	
	  Channel
		.fromPath(params.input)
		.map { path -> tuple(getSampleID(path.getName()), path) }
		.set {reads_ch}
	} else{
	
		Channel
		  .fromFilePairs(params.input)
		  .set {reads_ch}
	}
}

notSingleSample = !params.singleSample



if(params.singleEnd){
	process trimSE {
		tag "Trimmomatic on ${pair_id}"
		label 'med_mem'

		input:
		tuple val(pair_id), path(reads) from reads_ch
		
		output:
		path("${pair_id}_trim.log") into trimmomaticLogs_ch
		tuple pair_id, path("${pair_id}*.fastq.gz") into trimmedReads_ch, tReadsFqc_ch

		//TODO: add -threads $task.cpus
		//TODO: look into fastp as alternative trimmer to handle poly-G artefacts
		script:
		"""
		trimmomatic SE \
		  -threads $task.cpus \
		  ${reads[0]} \
		  ${pair_id}_trim_1P \
		  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 2> ${pair_id}_trim.log
		
		mv ${pair_id}_trim_1P ${pair_id}_trim_R1.fastq
		gzip ${pair_id}_trim_R1.fastq
		"""
	}

	process fastqcSE {
		
		tag "FASTQC on ${sample_id}"
		label 'small_mem'
		
		input:
		tuple val(sample_id), path(reads) from tReadsFqc_ch

		output:
		path("fastqc_${sample_id}_logs") into fastqc_ch

		script:
		"""
		mkdir fastqc_${sample_id}_logs
		fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
		"""  
	} 

	process bowtieAlignSE {
		tag "Aliging $pair_id to ${params.bt2_index}"
		label 'big_mem'

		input:
		val(idx) from params.bt2_index
		tuple val(pair_id), path(reads) from trimmedReads_ch

		output:
		path("${pair_id}_bt2.log") into bt2Logs_ch
		tuple pair_id, file("${pair_id}_iSort.bam"), file("${pair_id}_iSort.bam.bai") into bt2Bam_ch

		script:
		"""
		bowtie2 -p $task.cpus -x ${idx} --no-mixed --no-unal --no-discordant --local --very-sensitive-local -X 1000 -k 4 --mm -U ${reads} 2> ${pair_id}_bt2.log | samtools view -bS -q 30 - > ${pair_id}_init.bam
		samtools sort -@ $task.cpus ${pair_id}_init.bam > ${pair_id}_iSort.bam
		samtools index ${pair_id}_iSort.bam
		"""

	}


	process rmDupSE {
	   tag "Removing Dupes ${sampleID}"
	   label 'med_mem'
	   
	   input:
	   tuple val(sampleID), path(bam), path(index) from bt2Bam_ch
	   
	   output:
	   path("${sampleID}_dups.log") into picardDupStats_ch
	   tuple sampleID, file("${sampleID}_rmDup.bam") into rmDupBam_ch, idxStats_ch
	   
	   script:
	   """
	   picard MarkDuplicates VERBOSITY=WARNING \
		INPUT=${bam} OUTPUT=${sampleID}_rmDup.bam \
		METRICS_FILE=${sampleID}_dups.log \
		REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
	   """
	}

	process getIDXStatsSE {
	   tag "get mapping stats"
	   label 'small_mem'
	   
	   input:
	   tuple val(sampleID), path(bam) from idxStats_ch
	   
	   output:
	   path("${sampleID}_idxStats.log") into idxLog_ch
	   
	   script:
	   """
	   sambamba index ${bam}
	   samtools idxstats ${bam} > ${sampleID}_idxStats.log
	   """
	}

	process finalFilterSE {
	   tag "Removing chrM and BL"
	   label 'big_mem'
	   publishDir "$params.outdir/finalBam", mode: 'copy', pattern: "${sampleID}_final.bam"
	   
	   input:
	   path(blacklist) from params.blacklist
	   tuple val(sampleID), path(bam) from rmDupBam_ch
	   
	   
	   output:
	   tuple val(sampleID), file("${sampleID}_final.bam") finalBam_ch
	   file("${sampleID}_final.bam") into forPCA_ch, forBEPImage_ch
	   val(sampleID) into names_ch
	   
	   script:
	   """
	   samtools index ${bam}
	   export CHROMOSOMES=\$(samtools view -H ${bam} | grep '^@SQ' | cut -f 2 | grep -v -e _ -e chrM -e 'VN:' | sed 's/SN://' | xargs echo)
	   samtools view -b -h -f 3 -F 4 -F 256 -F 1024 -F 2048 -q 30 ${bam} \$CHROMOSOMES > tmp.bam
	   bedtools subtract -A -a tmp.bam -b ${blacklist} | samtools sort -@ $task.cpus - > ${sampleID}_final.bam
	   """
	}

	process makeBigwigSE{

		tag "Creating ${sampleID} bigwig"
		publishDir "$params.outdir/bigwig", mode: 'copy'
		label 'big_mem'

		input:
		tuple val(sampleID), file(finalBam) from finalBam_ch
		
		output:
		tuple val(sampleID), file("${sampleID}_CPMnorm.bw") into bigwig_ch, bigwig2_ch, bigwig3_ch
		val(sampleID) into labels_ch
		file("${sampleID}_CPMnorm.bw") into forGEnrichPlot_ch

		//TODO: add -p $task.cpus
		script:
		"""
		sambamba index $finalBam
		bamCoverage -p $task.cpus --bam ${finalBam} -o ${sampleID}_CPMnorm.bw -bs 10 --smoothLength 50 --normalizeUsing CPM --ignoreForNormalization chrX chrY  --skipNonCoveredRegions 
		"""
	}
	
	if (notSingleSample) {
		process  plotPCASE {
			tag "Creating bin-based Multi-Bam Summary"
			publishDir "$params.outdir/results", mode: 'copy'
		label 'massive_mem'
		
			input:
			path(files) from forPCA_ch.collect()
			val(name) from params.name
			val(pcaTitle) from params.PCATitle

			output:
			file("${name}_PCA.png") into results_ch

			//TODO: add -p $task.cpus
			script:
			"""
			for i in ${files}
			do
			  sambamba index \$i
			done
			multiBamSummary bins -p $task.cpus -b $files -o ${name}_matrix.npz --smartLabels 
			plotPCA \
			  -in ${name}_matrix.npz \
			  -o ${name}_PCA.png \
			  -T "$pcaTitle"
			"""

		}
	}




	process multiqcSE {
		publishDir "$params.outdir/results", mode:'copy'
		label 'small_mem'

		input:
		path('*') from fastqc_ch
		  .mix(idxLog_ch)
		  .mix(picardDupStats_ch)
		  .mix(trimmomaticLogs_ch)
		  .mix(bt2Logs_ch)
		  .collect()
		
		output:
		path('multiqc_report.html')

		script:
		"""
		multiqc .
		"""
	}
} else {
	process trim {
		tag "Trimmomatic on ${pair_id}"
		label 'med_mem'

		input:
		tuple val(pair_id), path(reads) from reads_ch
		
		output:
		path("${pair_id}_trim.log") into trimmomaticLogs_ch
		tuple pair_id, path("${pair_id}*.fastq.gz") into trimmedReads_ch, tReadsFqc_ch

		//TODO: add -threads $task.cpus
		//TODO: look into fastp as alternative trimmer to handle poly-G artefacts
		script:
		"""
		trimmomatic PE \
		  -threads $task.cpus \
		  ${reads[0]} \
		  ${reads[1]} \
		  -baseout ${pair_id}_trim \
		  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 2> ${pair_id}_trim.log
		
		mv ${pair_id}_trim_1P ${pair_id}_trim_R1.fastq
		mv ${pair_id}_trim_2P ${pair_id}_trim_R2.fastq
		gzip ${pair_id}_trim_R1.fastq
		gzip ${pair_id}_trim_R2.fastq
		"""
	}

	process fastqc {
		
		tag "FASTQC on ${sample_id}"
		label 'small_mem'
		
		input:
		tuple val(sample_id), path(reads) from tReadsFqc_ch

		output:
		path("fastqc_${sample_id}_logs") into fastqc_ch

		script:
		"""
		mkdir fastqc_${sample_id}_logs
		fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
		"""  
	} 

	process bowtieAlign {
		tag "Aliging $pair_id to ${params.bt2_index}"
		label 'big_mem'

		input:
		val(idx) from params.bt2_index
		tuple val(pair_id), path(reads) from trimmedReads_ch

		output:
		path("${pair_id}_bt2.log") into bt2Logs_ch
		tuple pair_id, file("${pair_id}_iSort.bam"), file("${pair_id}_iSort.bam.bai") into bt2Bam_ch

		//TODO: add -p $task.cpus
		script:
		"""
		bowtie2 -p $task.cpus -x ${idx} --no-mixed --no-unal --no-discordant --local --very-sensitive-local -X 1000 -k 4 --mm -1 ${reads[0]} -2 ${reads[1]} 2> ${pair_id}_bt2.log | samtools view -bS -q 30 - > ${pair_id}_init.bam
		
		samtools sort -@ $task.cpus ${pair_id}_init.bam > ${pair_id}_iSort.bam
		samtools index ${pair_id}_iSort.bam
		"""

	}


	process rmDup {
	   tag "Removing Dupes ${sampleID}"
	   label 'med_mem'
	   
	   input:
	   tuple val(sampleID), path(bam), path(index) from bt2Bam_ch
	   
	   output:
	   path("${sampleID}_dups.log") into picardDupStats_ch
	   tuple sampleID, file("${sampleID}_rmDup.bam") into rmDupBam_ch, idxStats_ch
	   
	   script:
	   """
	   picard MarkDuplicates VERBOSITY=WARNING \
		INPUT=${bam} OUTPUT=${sampleID}_rmDup.bam \
		METRICS_FILE=${sampleID}_dups.log \
		REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
	   """
	}

	process getIDXStats {
	   tag "get mapping stats"
	   label 'small_mem'
	   
	   input:
	   tuple val(sampleID), path(bam) from idxStats_ch
	   
	   output:
	   path("${sampleID}_idxStats.log") into idxLog_ch
	   
	   script:
	   """
	   sambamba index ${bam}
	   samtools idxstats ${bam} > ${sampleID}_idxStats.log
	   """
	}

	process finalFilterPE {
	   tag "Removing chrM and BL"
	   label 'big_mem'
	   publishDir "$params.outdir/finalBam", mode: 'copy', pattern: "${sampleID}_final.bam"
	   
	   input:
	   path(blacklist) from params.blacklist
	   tuple val(sampleID), path(bam) from rmDupBam_ch
	   
	   
	   output:
	   tuple val(sampleID), file("${sampleID}_final.bam") into CISM_ch, finalBam_ch
	   file("${sampleID}_final.bam") into forPCA_ch, forBEPImage_ch
	   val(sampleID) into names_ch
	   
	   script:
	   """
	   samtools index ${bam}
	   export CHROMOSOMES=\$(samtools view -H ${bam} | grep '^@SQ' | cut -f 2 | grep -v -e _ -e chrM -e 'VN:' | sed 's/SN://' | xargs echo)
	   samtools view -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 ${bam} \$CHROMOSOMES > tmp.bam
	   bedtools subtract -A -a tmp.bam -b ${blacklist} | samtools sort -@ $task.cpus - > ${sampleID}_final.bam
	   """
	}

	process collectInserts {
	   tag "generating Insert Sizes"
	   label 'big_mem'
	   
	   input:
	   tuple val(sampleID), path(bam) from CISM_ch
	   
	   output:
	   path("${sampleID}_insertSizes.log") into picardISStats_ch
	   
	   script:
	   """
	   samtools index ${bam}
	   picard CollectInsertSizeMetrics \
		I=${bam} \
		O=${sampleID}_insertSizes.log \
		H=${sampleID}_insertHist.pdf 
	   """
	}

	process makeBigwig{

		tag "Creating ${sampleID} bigwig"
		publishDir "$params.outdir/bigwig", mode: 'copy'
		label 'big_mem'

		input:
		tuple val(sampleID), file(finalBam) from finalBam_ch
		
		output:
		tuple val(sampleID), file("${sampleID}_CPMnorm.bw") into bigwig_ch, bigwig2_ch, bigwig3_ch
		val(sampleID) into labels_ch
		file("${sampleID}_CPMnorm.bw") into forGEnrichPlot_ch

		//TODO: add -p $task.cpus
		script:
		"""
		sambamba index $finalBam
		bamCoverage -p $task.cpus --bam ${finalBam} -o ${sampleID}_CPMnorm.bw -bs 10 --extendReads --smoothLength 50 --normalizeUsing CPM --ignoreForNormalization chrX chrY  --skipNonCoveredRegions 
		"""
	}

	BEFPDF_ch = names_ch.toSortedList()
	sortedNamedBam = forBEPImage_ch.toSortedList()


	process generateGlobalFragmentPDF {
		tag "Creating Summary Fragment Histograms"
		publishDir "$params.outdir/results", mode: 'copy'
		label 'med_mem'

		input:
		path(files) from sortedNamedBam
		val(labels) from BEFPDF_ch
		val(name) from params.name

		output:
		file( "${name}_PEFragHist-all.pdf" )
		
		//TODO: add -p $task.cpus
		script:
		"""
		for i in $files; do
		  sambamba index \$i
		done
		bamPEFragmentSize -p $task.cpus -b ${files} -o ${name}_PEFragHist-all.pdf --samplesLabel ${labels}
		"""
	}
	
	if (notSingleSample) {
    process  plotPCA {
        tag "Creating bin-based Multi-Bam Summary"
        publishDir "$params.outdir/results", mode: 'copy'
	label 'massive_mem'
    
        input:
        path(files) from forPCA_ch.collect()
        val(name) from params.name
        val(pcaTitle) from params.PCATitle

        output:
        file("${name}_PCA.png") into results_ch

		//TODO: add -p $task.cpus
        script:
        """
        for i in ${files}
        do
          sambamba index \$i
        done
        multiBamSummary bins -p $task.cpus -b $files -o ${name}_matrix.npz --smartLabels --extendReads
        plotPCA \
          -in ${name}_matrix.npz \
          -o ${name}_PCA.png \
          -T "$pcaTitle"
        """
		}
	}

	process multiqc {
		publishDir "$params.outdir/results", mode:'copy'
		label 'small_mem'

		input:
		path('*') from fastqc_ch
		  .mix(idxLog_ch)
		  .mix(picardISStats_ch)
		  .mix(picardDupStats_ch)
		  .mix(trimmomaticLogs_ch)
		  .mix(bt2Logs_ch)
		  .collect()
		
		output:
		path('multiqc_report.html')

		script:
		"""
		multiqc .
		"""
	}
}







process computeMatrixDefault {
    tag "${sampleID} generating gene-wide TSS and GB profile matrices"
    label 'med_mem'
    
    input:
    tuple val(sampleID), file(bigwig) from bigwig_ch
    path(genes) from params.genesList

    output:
    tuple val(sampleID), path("${sampleID}_rpMat.npz") into tssMatrixGW_ch
    tuple val(sampleID), path("${sampleID}_srMat.npz") into profileMatrixGW_ch
    file("${sampleID}_rpMat.npz") into tssMatrixglob_ch
    file("${sampleID}_srMat.npz") into srMatrixglob_ch

	//TODO: add -p $task.cpus
    script:
    """
    computeMatrix reference-point -p $task.cpus  -S $bigwig -R $genes -o ${sampleID}_rpMat.npz -b 3000 -a 3000 --missingDataAsZero --samplesLabel ${sampleID}
    computeMatrix scale-regions -p $task.cpus -S $bigwig -R $genes -o ${sampleID}_srMat.npz -m 8000 -b 3000 -a 3000 --missingDataAsZero --samplesLabel $sampleID
    """
}




process generateEnrichPlots {
    tag "${sampleID} TSS and Gene-body Enrichment"
    publishDir "${params.outdir}/results/${sampleID}", mode: 'copy', pattern: "*.pdf"
    label 'small_mem'

    input:
    tuple val(sampleID), file(matrix) from tssMatrixGW_ch
    tuple val(sampleID2), file(matrix2) from profileMatrixGW_ch

    output: 
    file("${sampleID}_tssEnrich.pdf") 
    file("${sampleID}_geneBodyProfile.pdf")

    script:
    """
    plotHeatmap -m $matrix -o ${sampleID}_tssEnrich.pdf
    plotHeatmap -m $matrix2 -o ${sampleID}_geneBodyProfile.pdf
    """
}


process makeGlobalEnrichPlots {
    tag "Project: ${name} TSS and Gene Body Plots"
    publishDir "$params.outdir/results", mode: 'copy', pattern: "*.pdf"
    label 'largeStore'
    
    input:
    val(name) from params.name
    path(files) from tssMatrixglob_ch.toSortedList().collect()
    path(files2) from srMatrixglob_ch.toSortedList().collect()

    output:
    file("${name}_TSSPlot.pdf")
    file("${name}_GeneBodyPlot.pdf")
    
    script:
    """
    
    computeMatrixOperations cbind -m $files -o ${name}_TSSmat.npz
    plotHeatmap -m ${name}_TSSmat.npz -o ${name}_TSSPlot.pdf
    
    computeMatrixOperations cbind -m $files2 -o ${name}_GeneBodymat.npz
    plotHeatmap -m ${name}_GeneBodymat.npz -o ${name}_GeneBodyPlot.pdf
    """


}






if(params.addBEDFilesProfile) {
    Channel
      .fromPath(params.addBEDFilesProfile)
	  .splitCsv(header:false)
	  .map{ row -> tuple(row[0], file(row[1])) }
	  .set { extraBEDs_ch }

    extraBEDs_ch
	  .combine(bigwig2_ch)
	  .set { totalExtraBed_ch }

   
    process computeMatExtra {
        tag "Compute Matrix for ${sampleID} on extra BED file: ${extraBEDName}"
        label 'med_mem'
        
        input:
        tuple val(extraBEDName), path(BED), val(sampleID), path(bigwig) from totalExtraBed_ch

        output:
        tuple val(extraBEDName), val(sampleID), file("${sampleID}-${extraBEDName}_profile.npz") into addBEDMatTuple_ch
        tuple val(extraBEDName), file("${sampleID}-${extraBEDName}_profile.npz") into addBEDMatTupleGlobal_ch
    
		//TODO: add -p $task.cpus
        script:
        """
        computeMatrix scale-regions -p $task.cpus -S $bigwig -R $BED -b 3000 -m 8000 -a 3000 --missingDataAsZero --samplesLabel ${sampleID} -o ${sampleID}-${extraBEDName}_profile.npz
        """
    }
    


    process generateExtraBEDProfiles {
        tag "Visualizing read density for ${rName} on sample ${sName}"
        publishDir "$params.outdir/results/extraBED/${sName}", mode: 'copy'
	    label 'small_mem'

        input:
        tuple val(rName), val(sName), path(mat) from addBEDMatTuple_ch
        output:
        file("${sName}-${rName}_profile.pdf")
 
        script:
        """
        plotHeatmap -m $mat -z $rName -o "${sName}-${rName}_profile.pdf"
        """

    } 

    addBEDMatTupleGlobal_ch
      .groupTuple()
      .set { mixedExtraBEDsGT_ch }

    process generateGlobalExtraBED {
        tag "Combining profile plots for ${rName}"
        publishDir "$params.outdir/results/extraBED", mode: 'copy', pattern: "*.pdf"
	    label 'largeStore'
    
        input:
        tuple val(rName), path(mats) from mixedExtraBEDsGT_ch
        
        output:
        file("${rName}_profileAll.pdf")

        script:
        """
        computeMatrixOperations cbind -m ${mats} -o ${rName}_gMat.npz
        plotHeatmap -m ${rName}_gMat.npz -z ${rName} -o ${rName}_profileAll.pdf
        """
    }

}


if(params.addBEDFilesRefPoint) {
    Channel
      .fromPath(params.addBEDFilesRefPoint)
	  .splitCsv(header:false)
	  .map{ row -> tuple(row[0], file(row[1]), row[2], row[3]) }
	  .set { extraBEDs2_ch }

    extraBEDs2_ch
	  .combine(bigwig3_ch)
	  .set { totalExtraBed2_ch }
	  
	

   
    process computeMatExtraRP {
        tag "Compute Matrix for ${sampleID} on extra BED file: ${extraBEDName}"
	    label 'med_mem'
        
        input:
        tuple val(extraBEDName), path(BED), val(range), val(pointLabel), val(sampleID), path(bigwig) from totalExtraBed2_ch

        output:
        tuple val(extraBEDName), val(sampleID), file("${sampleID}-${extraBEDName}_refPoint.npz") into addBEDMatTuple2_ch
        tuple val(extraBEDName), file("${sampleID}-${extraBEDName}_refPoint.npz") into addBEDMatTupleGlobal2_ch
    
		//TODO: add -p $task.cpus
        script:
        """
        computeMatrix reference-point -p $task.cpus -S $bigwig -R $BED -b ${range} -a ${range} --missingDataAsZero --samplesLabel ${sampleID} -o ${sampleID}-${extraBEDName}_refPoint.npz
        """
    }
    


    process generateExtraBEDRP {
        tag "Visualizing read density for ${rName} on sample ${sName}"
        publishDir "$params.outdir/results/extraBED/${sName}", mode: 'copy'
	    label 'small_mem'

        input:
        tuple val(rName), val(sName), path(mat) from addBEDMatTuple2_ch
        output:
        file("${sName}-${rName}_refPoint.pdf")
 
        script:
        """
        plotHeatmap -m $mat --refPointLabel $rName -o "${sName}-${rName}_refPoint.pdf"
        """

    } 

    addBEDMatTupleGlobal2_ch
      .groupTuple()
      .set { mixedExtraBEDsGT2_ch }

    process generateGlobalExtraBEDRP {
        tag "Combining profile plots for ${rName}"
        publishDir "$params.outdir/results/extraBED", mode: 'copy', pattern: "*.pdf"
	    label 'largeStore'
    
        input:
        tuple val(rName), path(mats) from mixedExtraBEDsGT2_ch
        
        output:
        file("${rName}_refPointAll.pdf")

        script:
        """
        computeMatrixOperations cbind -m ${mats} -o ${rName}_gMat.npz
        plotHeatmap -m ${rName}_gMat.npz -z ${rName} -o ${rName}_refPointAll.pdf
        """
    }

}    
