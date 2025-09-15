# singlecell-analysis-pipeline

Welcome to our collaborative pipeline project!
This repository is created and maintained by Junbeom and Sohyeong as part of our effort to learn and build a bioinformatics workflow.
Our goal is to write and refine scripts step by step, and then connect them into a final automated and interactive single cell ATAC-seq pipeline. The pipeline will take SRR raw data, convert it into the proper format for Cell Ranger, run the complete Cell Ranger processing, and finally use the processed results for peak-level analysis.
Through this project, we not only aim to build a useful tool for single-cell data analysis, but also to gain hands-on experience in collaborative coding with GitHub, version control, and teamwork.

---

## Folder structure

```
singlecell-analysis-pipeline/
├── config/      # Sample sheet and parameter files
├── scripts/     # Step-wise scripts (sh, R, Python) -> 여기에 개별 코드를 업로드
├── results/     # Final outputs (ignored in Git, use .gitkeep only)
├── logs/        # Log files from each step (ignored in Git)
├── work/        # Intermediate files, cache (ignored in Git)
├── tests/       # Small test data and test scripts
├── README.md    # Project overview
├── LICENSE      # License info (MIT)
└── .gitignore   # Ignore rules for data, results, cache
└── pipeline.sh  # 최종 파이프라인 실행 코드
```

---

## How to upload your code to GitHub

> This is a simple guide for beginners who have never used GitHub before.

1. **Clone the repository (first time only)**
   - Go to your wanted path in your server or local path.
   - Type as follows.
   ```bash
   git clone https://github.com/seq-jchoi-bio/singlecell-analysis-pipeline.git
   cd singlecell-analysis-pipeline
   ```

2. **Check the folder structure**
   ```bash
   ls
   # You should see: config/ scripts/ results/ logs/ work/ tests/ ...
   ```

3. **Create your NEW file!**
   - Put your scripts into `scripts/`
   - Example: `scripts/test.sh`
   ```bash
   nano scripts/test.sh
   # (write your code here, save & exit)
   ```

4. **Stage → Commit → Push**
   ```bash
   git status
   git add scripts/test.sh
   git commit -m "Add test script"
   git push origin main
   ```

   - Go to the repository page in your browser  
   - Confirm that your file appears in the correct folder

5. **Modify your script!**
   - Example: `scripts/test.sh`
   ```bash
   nano scripts/test.sh
   # (modify your code here, save & exit)
   ```

6. **Stage → Commit → Push**
   ```bash
   git status
   git add scripts/test.sh
   git commit -m "Update test code: fix file"
   git push origin main
   ```

   - Go to the repository page in your browser  
   - Confirm that your file appears in the correct folder

## Tips
- Always run `git status` before commit to see what changed.
- Use clear commit messages:
  - "Add new QC script" (for new file)  
  - "Update peak-call script: change MACS2 options" (for modification)
