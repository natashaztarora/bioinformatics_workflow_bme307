# Installation Instructions

All students must bring their own laptop and install the required software **before the start of the practicals**.

We will provide a **troubleshooting session on Tuesday (15.09.2026)** after the lecture for those who encounter installation issues.

## Required Software

1. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) — to inspect the quality of raw sequencing reads.
2. [R](https://rstudio-education.github.io/hopr/starting.html) — the programming language we will use for downstream microbiome analysis.
3. [R Studio](https://posit.co/products/open-source/rstudio) — the interface we will use to work with R.

If you encounter installation problems, please contact the course instructors before the practical at **[bioinformatics.bme307@gmail.com](mailto:bioinformatics.bme307@gmail.com)**.

---

## 💻 Windows Users 

### 1. FastQC

1. Go to the [FastQC download page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. Download the **Windows/Linux zip file** for the current FastQC release.
3. Extract the downloaded `.zip` file.
4. Open the extracted FastQC folder.
5. Start FastQC using the Windows launcher included in the folder.

!!! note
    FastQC requires Java. If FastQC does not start correctly, please contact the course instructors.

### 2. R 

1. Go to [CRAN](https://cran.r-project.org/).
2. Click **Download R for Windows**.
3. Click **base** and then download the current R installer.
4. Open the downloaded installer and keep the default installation options.

### 3. RStudio Desktop

1. Go to the [RStudio Desktop download page](https://posit.co/download/rstudio-desktop/).
2. Download the **free RStudio Desktop** installer for Windows.
3. Run the installer using the default options.
4. Open RStudio.

To check that R is working, type the following in the **Console** and press Enter:

```r
1 + 1
```

You should see:

```text
[1] 2
```

---

## 🍎 Mac Users 

### 1. FastQC

1. Go to the [FastQC download page](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. Download the **Mac DMG image** for the current FastQC release.
3. Open the downloaded `.dmg` file and install FastQC.
4. Open FastQC once to check that it launches correctly.

!!! note
    If macOS blocks the application the first time you open it, go to **System Settings → Privacy & Security** and allow the application to open.

### 2. R

1. Go to [CRAN](https://cran.r-project.org/).
2. Click **Download R for macOS**.
3. Download the installer that matches your Mac:
    - **Apple Silicon (`arm64`)** for M1/M2/M3/M4/M5 Macs.
    - **Intel (`x86_64`)** for older Intel Macs.
4. Open the `.pkg` file and follow the installation instructions.

If you are unsure which Mac you have, click **Apple menu → About This Mac** and look at **Chip** or **Processor**.

### 3. RStudio Desktop

1. Go to the [RStudio Desktop download page](https://posit.co/download/rstudio-desktop/).
2. Download the **free RStudio Desktop** version for macOS.
3. Open the `.dmg` file and drag RStudio into **Applications**.
4. Open RStudio.

To check that R is working, type the following in the **Console** and press Enter:

```r
1 + 1
```

You should see:

```text
[1] 2
```

---

## 🐧 Linux Users

The instructions below are for **Ubuntu/Debian-based distributions**. If you use another Linux distribution, please use the corresponding packages for your system.

### 1. FastQC

Open a terminal and run:

```bash
sudo apt update
sudo apt install fastqc
```

Check the installation with:

```bash
fastqc --version
```

### 2. R

Install R with:

```bash
sudo apt update
sudo apt install r-base
```

Check the installation with:

```bash
R --version
```

### 3. RStudio Desktop

1. Go to the [RStudio Desktop download page](https://posit.co/download/rstudio-desktop/).
2. Download the **free RStudio Desktop** package corresponding to your Linux distribution.
3. Install the downloaded package following the instructions on the download page.
4. Open RStudio.

In the RStudio **Console**, run:

```r
1 + 1
```

You should see:

```text
[1] 2
```

---

## Optional: QIIME 2

??? info "QIIME 2 installation (optional — not required for the course)"
    **QIIME 2 is not required for the practical.**

    During the course, the QIIME 2 processing steps will be explained, but the processing itself will already have been run for you. You will receive the files needed to continue the analysis in RStudio.

    If you would like to reproduce the complete workflow independently after the course, follow the official installation instructions for the current QIIME 2 Amplicon distribution:

    [QIIME 2 installation guide](https://amplicon-docs.qiime2.org/en/latest/how-to-guides/install.html)

---

## Before the practical: quick checklist

Please make sure that:

- [ ] FastQC opens successfully.
- [ ] R is installed.
- [ ] RStudio Desktop opens successfully.
- [ ] Running `1 + 1` in the RStudio Console returns `[1] 2`.
- [ ] You understand that **QIIME 2 does not need to be installed for the course**.
