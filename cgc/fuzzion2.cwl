class: CommandLineTool
cwlVersion: v1.2
baseCommand: []
inputs:
  - id: pattern
    type: File
    inputBinding:
      prefix: '-pattern='
      separate: false
      shellQuote: false
      position: 0
    'sbg:fileTypes': .txt
  - id: rank
    type: File
    inputBinding:
      prefix: '-rank='
      separate: false
      shellQuote: false
      position: 1
    'sbg:fileTypes': .krt
  - id: fastq1
    type: File
    inputBinding:
      prefix: '-fastq1='
      separate: false
      shellQuote: false
      position: 2
    'sbg:fileTypes': '.fq, .fastq, .fq.gz, .fastq.gz'
  - id: fastq2
    type: File?
    inputBinding:
      prefix: '-fastq2='
      separate: false
      shellQuote: false
      position: 3
    'sbg:fileTypes': '.fq, .fastq, .fq.gz, .fastq.gz'
outputs:
  - id: example_out
    type: stdout
    outputBinding:
      glob: output.txt
label: fuzzion2
requirements:
  - class: ShellCommandRequirement
  - class: ResourceRequirement
    ramMin: 5000
  - class: DockerRequirement
    dockerPull: 'cgc-images.sbgenomics.com/cjrash/stjude/fuzzion2:latest'
stdout: output.txt
