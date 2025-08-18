# Installation Instructions

All students must bring their own laptop and install the required software **before the start of the practicals**.

We will provide a **troubleshooting session on Tuesday (16.9.2025)** after the lecture for those who encounter installation issues.

## Required Software

1. [Docker](https://www.docker.com/)
2. [QIIME2](https://qiime2.org/)
3. [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
4. [R + RStudio](https://rstudio-education.github.io/hopr/starting.html)

If you have questions, please email us at **[bioinformatics.bme307@gmail.com](mailto:bioinformatics.bme307@gmail.com)**.

---

## 💻 Windows Users (recommended: Docker-based QIIME2)

### 1. Docker + QIIME2

!!! warning "Check if Virtualization is enabled"
    - Press **Ctrl + Alt + Del → Task Manager → Performance → CPU**
    - Check if **Virtualization: Enabled**
    - If not, follow [these BIOS instructions](https://wiki.2n.com/faqac/en/virtualizace-vt-x-amd-v-povoleni-virtualizace-na-vasem-pocitaci-pro-spusteni-2n-access-commander-100572533.html) (steps differ by brand: Acer, Asus, Dell, HP, Lenovo, Sony, Toshiba).

1. Install Docker Desktop for Windows: [Download here](https://docs.docker.com/desktop/install/windows-install/)

   * Run the `.exe` file and select **Enable WSL 2 Features** when prompted.

2. If you see:
   **“Docker Desktop requires Windows …”**

   * Update Windows [here](https://www.microsoft.com/en-us/software-download/windows10) and reinstall Docker.

3. If you encounter **WSL2 errors**, install the [Linux kernel update package](https://wslstorestorage.blob.core.windows.net/wslblob/wsl_update_x64.msi).

4. Start Docker Desktop.

   * The whale 🐳 icon in the status bar should stay steady.

5. Test Docker with:

   ```bash
   docker run hello-world
   ```

6. Download QIIME2:

   ```bash
   docker pull quay.io/qiime2/amplicon:2025.7
   ```

7. Test QIIME2:

   ```bash
   docker run -v ${PWD}:/data -it quay.io/qiime2/amplicon:2025.7 qiime info
   ```

---

### 2. FastQC

1. Download FastQC for Windows: [FastQC v0.12.1](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. Extract the `.zip` file.
3. Run `fastqc.exe`.

   * If Java errors occur, install [Java JDK 21](https://www.oracle.com/java/technologies/downloads/#jdk21-windows).

---

### 3. R + RStudio

1. Install [R for Windows](http://cran.r-project.org/bin/windows/base/release.htm).
2. Install [RStudio for Windows](https://posit.co/download/rstudio-desktop/).
3. Open RStudio to confirm installation works.

---

## 🍎 Mac Users (Intel & Apple Silicon M1–M4)

### 1. Install Conda

* Recommended: [Miniconda](https://docs.anaconda.com/miniconda/install/#macos-terminal-installer)
* Open Terminal and verify:

  ```bash
  conda --version
  ```

---

### 2. QIIME2

#### Intel Macs (x86\_64)

```bash
conda update conda
conda env create --name qiime2-amplicon-2025.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
conda activate qiime2-amplicon-2025.7
```

---

#### Apple Silicon Macs (M1–M4, arm64)

Use the `osx-64` builds via Rosetta:

```bash
CONDA_SUBDIR=osx-64 conda env create \
  --name qiime2-amplicon-2025.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-macos-latest-conda.yml

conda activate qiime2-amplicon-2025.7
conda config --env --set subdir osx-64
```

---

#### Test installation

```bash
qiime info
```

If successful, this will print information about the QIIME2 version.

---

### 3. FastQC

1. Download the `.dmg`: [FastQC v0.12.1](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.dmg).
2. If macOS blocks the app, install [Java JDK 21 for Mac](https://www.oracle.com/java/technologies/downloads/#jdk21-mac).

If FastQC still doesn’t launch:

* Download the [source code](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip).
* Unzip → open Terminal → navigate to the unzipped folder → run:

```bash
./fastqc
```

---

#### 🔗 Make FastQC executable from anywhere

After unzipping, create a **symbolic link** so `fastqc` works globally.

**Intel Macs (x86\_64):**

```bash
sudo ln -s ~/Downloads/FastQC/fastqc /usr/local/bin/fastqc
sudo chmod +x /usr/local/bin/fastqc
```

**Apple Silicon Macs (M1–M4, arm64):**

```bash
sudo ln -s ~/Downloads/FastQC/fastqc /opt/homebrew/bin/fastqc
sudo chmod +x /opt/homebrew/bin/fastqc
```

*(adjust the path if your FastQC folder is elsewhere)*

---

#### Test the installation

Run:

```bash
fastqc -h
```

If installed correctly, this will print the FastQC help menu.

---


---

### 4. R + RStudio

1. Install [R for macOS](https://cran.r-project.org/bin/macosx/).
2. Install [RStudio for macOS](https://posit.co/download/rstudio-desktop/).
3. Open RStudio to confirm installation.

---

## 🐧 Linux Users

### 1. Docker (recommended for simplicity)

Install Docker:

```bash
sudo apt-get update
sudo apt-get install docker.io
```

Test Docker:

```bash
docker run hello-world
```

Download and test QIIME2:

```bash
docker pull quay.io/qiime2/amplicon:2025.7
docker run -v ${PWD}:/data -it quay.io/qiime2/amplicon:2025.7 qiime info
```

---

### 2. FastQC

```bash
sudo apt-get install fastqc
```

Or download from [FastQC site](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

---

### 3. R + RStudio

```bash
sudo apt-get install r-base
```

Download RStudio: [Linux RStudio](https://posit.co/download/rstudio-desktop/).
