# Table of Contents

| Section | Page |
|---|---:|
| Abstract | 4 |
| Acknowledgements | 5 |
| **1. Introduction** | 6 |
| &nbsp;&nbsp;&nbsp;&nbsp;1.1 IDPs and coupled folding-and-binding | 6 |
| &nbsp;&nbsp;&nbsp;&nbsp;1.2 Two candidate information sources for d2o prediction | 6 |
| &nbsp;&nbsp;&nbsp;&nbsp;1.3 Study design and contributions | 7 |
| **2. Methods** | 8 |
| &nbsp;&nbsp;&nbsp;&nbsp;2.1 Dataset and label | 8 |
| &nbsp;&nbsp;&nbsp;&nbsp;2.2 Features | 8 |
| &nbsp;&nbsp;&nbsp;&nbsp;2.3 Sequence-redundancy control and cross-validation | 8 |
| &nbsp;&nbsp;&nbsp;&nbsp;2.4 Factorial design and classifier | 9 |
| &nbsp;&nbsp;&nbsp;&nbsp;2.5 Metrics and statistical analysis | 9 |
| &nbsp;&nbsp;&nbsp;&nbsp;2.6 Extensions | 10 |
| &nbsp;&nbsp;&nbsp;&nbsp;2.7 Limitations | 11 |
| **3. Results** | 12 |
| &nbsp;&nbsp;&nbsp;&nbsp;3.1 Dataset and primary factorial | 12 |
| &nbsp;&nbsp;&nbsp;&nbsp;3.2 PLMs already contain GO Slim's discriminative content | 12 |
| &nbsp;&nbsp;&nbsp;&nbsp;3.3 GO context acts as a regularizer toward the class prior | 13 |
| &nbsp;&nbsp;&nbsp;&nbsp;3.4 Sequence representation is the strongest lever | 14 |
| &nbsp;&nbsp;&nbsp;&nbsp;3.5 Higher-resolution GO does not rescue the null | 15 |
| &nbsp;&nbsp;&nbsp;&nbsp;3.6 Robustness across classifier, PLM family, and region scale | 16 |
| **4. Discussion** | 18 |
| &nbsp;&nbsp;&nbsp;&nbsp;4.1 Principal finding | 18 |
| &nbsp;&nbsp;&nbsp;&nbsp;4.2 Why PLMs subsume GO Slim | 18 |
| &nbsp;&nbsp;&nbsp;&nbsp;4.3 Regularization toward the class prior | 19 |
| &nbsp;&nbsp;&nbsp;&nbsp;4.4 What the higher-resolution GO test told us | 19 |
| &nbsp;&nbsp;&nbsp;&nbsp;4.5 Positioning, limitations, future work | 19 |
| &nbsp;&nbsp;&nbsp;&nbsp;4.6 Biological implications | 20 |
| &nbsp;&nbsp;&nbsp;&nbsp;4.7 Conclusion | 21 |
| **References** | 23 |
| **Appendices** | 28 |
| &nbsp;&nbsp;&nbsp;&nbsp;A. Per-fold AUPRC for the disorder-pool factorial | 28 |
| &nbsp;&nbsp;&nbsp;&nbsp;B. ProstT5 replication | 28 |
| &nbsp;&nbsp;&nbsp;&nbsp;C. Region-level factorial | 28 |
| &nbsp;&nbsp;&nbsp;&nbsp;D. Higher-resolution GO factorial and strengthening pass | 29 |
| &nbsp;&nbsp;&nbsp;&nbsp;E. Reproducibility | 31 |
