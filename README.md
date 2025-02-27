**Management and rhizosphere microbial associations modulate genetic-driven nitrogen fate**
---

Mitra Ghotbi<sup>1,2*</sup>, Marjan Ghotbi<sup>3</sup>, Yakov Kuzyakov<sup>4,5</sup>, William R. Horwath<sup>6</sup>


**Affiliations:**

**Mitra Ghotbi**

- <sup>1</sup> Department of Natural Resources & Environmental Sciences, University of Illinois at Urbana-Champaign, IL, USA
- <sup>2</sup> Department of Biology, Middle Tennessee State University at Murfreesboro, TN, USA

**Marjan Ghotbi**

- <sup>3</sup> GEOMAR Helmholtz Centre for Ocean Research at Kiel, Germany

**Yakov Kuzyakov**

- <sup>4</sup> Department of Soil Science of Temperate Ecosystems and Department of Agricultural Soil Science, Georg-August-University of Göttingen at Göttingen, Germany
- <sup>5</sup> Peoples Friendship University of Russia (RUDN University), 117198 Moscow, Russia

**William R. Horwath**

- <sup>6</sup> Plant and Environmental Sciences Building, Department of Land, Air & Water Resources, University of California Davis, at Davis, CA, USA

**Corresponding Author:**

Mitra Ghotbi, PhD.

- Department of Biology, Middle Tennessee State University at Murfreesboro, TN, USA
- Email: mitra.ghotbi@mtsu.edu
- ORCID ID: [0000-0001-9185-9993](https://orcid.org/0000-0001-9185-9993)

---
![N Fate](https://github.com/mghotbi/Rhizosphere-Nitrogen-Fate/blob/Rhizosphere-nitrogen-fate/GA_AGEE%20(1).png)

**The R codes relevant to each figure have been published here to assist readers in navigating and reproducing them easily.**

# Prerequisites: Ensure that the required packages are installed and loaded.

```r
reqpkg <- c("phyloseq", "vegan", "microbiome", "DESeq2", 
            "plyr", "reshape2", "grid", "scales", "cluster", "ape", "dplyr","ggrepel",
            "igraph", "ggnet", "microbiomeutilities", "network", "SpiecEasi",
            "data.table", "decontam", "ggtext","microbiomeMarker", "devtools", "dada2", "ggplot2", "ggpubr",
            "agridat", "lme4", "rstatix", "emmeans", "lmerTest")

# Check against installed packages:
inpkg <- installed.packages()[, "Package"]
neededpkg <- reqpkg[!reqpkg %in% inpkg]

if (length(neededpkg) > 0) {
  cat("The following package(s) need to be installed:", paste(neededpkg, collapse = ", "), "\n")
  cat("Installing package(s)...\n")
  install.packages(neededpkg)
  cat("Package(s) installed successfully.\n")
} else {
  cat("All required packages are already installed.\n")
}

# Load required packages:
cat("Loading required packages...\n")
invisible(lapply(reqpkg, library, character.only = TRUE))
cat("All required packages loaded successfully.\n")
```

---

