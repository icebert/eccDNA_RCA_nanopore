## Base Calling, Mapping and QC of Nanopore sequencing

This snakemake based pipeline integrates:

1. base calling from raw fast5 files generated from Oxford Nanopore sequencer
2. adaptor trimming using porechop
3. reads mapping with minimap2
4. QC using [minionQC](https://github.com/roblanf/minion_qc), [NanoporeTools](https://github.com/Nucleomics-VIB/nanopore-tools) and [POREquality](https://github.com/carsweshau/POREquality)

The output mapped file is in [PAF](https://lh3.github.io/minimap2/minimap2.html#10) format.

#### Dependency

- [Python3](https://www.python.org)
- [R](https://www.r-project.org)
- [snakemake](https://snakemake.github.io)
- [guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis)
- [porechop](https://github.com/rrwick/Porechop)
- [minimap2](https://github.com/lh3/minimap2)

#### Run

To run the pipeline, first edit the config.json with actual file path, then run
using command:
```
snakemake -j 8 -p
```

