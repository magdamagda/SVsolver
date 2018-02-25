# SVsolver
#### Tool for detecting structural rearrangements in genom using long reads

### Usage
1. Prepare configuration file. Example of configuration file is in /example/config.ini
2. Run `./SVsolver/SVsolver.py <path to configuration file>`

### How to run example
1. Prepare coverage files for each contig using `samtools depth` and save them in example/coverage/
2. Create output directory example/output/
3. Run ./SVsolver/SVsolver.py /example/config.ini
