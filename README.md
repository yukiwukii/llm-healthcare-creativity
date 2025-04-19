# LLM Healthcare Creativity

## 🌏 Overview
This is the code repository for 'Benchmarking LLM Creativity in Healthcare and Drug Design' [(Full Paper)](Report.pdf).

## 💾 Installation and Preparation
1. First, clone this repository and `cd` into it.
```bash
# Clone the llm-healthcare-creativity repository
git clone https://github.com/yukiwukii/llm-healthcare-creativity.git

# Move into the repository directory
cd llm-healthcare-creativity
```

2. Install the necessary dependencies.
```bash
# Optional: Create a virtual environment
# conda create -n healthcreative python=3.11
# conda activate healthcreative

# Install the dependencies
pip install -r requirements.txt
```

3. Create a `.env` file consisting of your OpenAI API Key for GPT-as-a-judge.
```bash
OPENAI_API_KEY = 'your_api_key_here'
```

## 🚀 Running Experiments
### Medical QA Task
1. Run the LLM inference script.
```bash
ollama pull mistral, gemma2, llama3.1
python medical_qa_task/inference_for_medical_qa.py
```
2. Run GPT-as-a-judge script for both convergent and divergent creativity.
```bash
# For convergent creativity
python medical_qa_task/convergent_creativity_judge.py

# For divergent creativity
python medical_qa_task/divergent_creativity_judge.py
```
3. Analyse results using the given jupyter notebooks: `medical_qa_task/convergent_creativity_analysis` and `medical_qa_task/divergent_creativity_analysis`.

### Protein Generation Task
1. Run the LLM inference script.
```bash
ollama pull deepseek-r1:8b, gemma2, llama3.1
python protein_generation_task/inference_for_protein_generation.py
```
2. Run AlphaFold2 using [localcolabfold](https://github.com/YoshitakaMo/localcolabfold).
```bash
# Install localcolabfold
git clone https://github.com/YoshitakaMo/localcolabfold.git
cd localcolabfold
wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
bash install_colabbatch_linux.sh

# Add environment variable PATH
export PATH="/path/to/your/localcolabfold/colabfold-conda/bin:$PATH"

# Run AlphaFold2, inputdir is your generated fasta files from step 1
colabfold_batch inputdir outputdir/
```
3. Run PyRosetta Docking. We used ompT as the target protein, but you may use another protein, as long as you have its `.pdb` file. Move the generated `.pdb` files from step 2 into `pdb_files/protein`. 
```bash
python protein_generation_task/pyrosetta/directory_dock.py --receptor ompT_cleaned.pdb --ligand_dir pdb_files/protein --output_dir docking_results/protein
```
4. In the output directory, `docking_results/protein`, find the `docking_summary.tsv` file. In the same folder as `.pdb` files in step 2, find the `log.txt` file. Put them in the `./binding_analysis` folder. Run binding analysis.
```bash
python protein_generation_task/binding_analysis/binding_analysis.py
```
5. Run novelty analysis using BLAST.
```bash
python protein_generation_task/protein_novelty/blast.py
```
6. Analyse protein novelty results with HSSP.
```bash
python protein_generation_task/protein_novelty/hssp_analysis.py
```

## 💡 Samples
We provided results of the experiments in the paper at `./samples`. 
| **File Name**                                      | **Description**                                             |
| --------------------------------------------------- | -------------------------------------------------------- |
| `./convergent_creativity_score.csv`| [Medical QA Task] GPT-as-a-judge results for convergent creativity |
| `./divergent_creativity_score.csv`        | [Medical QA Task] GPT-as-a-judge results for divergent creativity            |
| `./generated_fasta_files`           |[Protein Generation Task] Generated FASTA files by the LLMs                           |
| `./pdb_files`| [Protein Generation Task] Structures of generated sequences after folding by AlphaFold2|
| `./docking_results` | [Protein Generation Task] Docking results of generated sequences after docking by PyRosetta |
| `./blast_results` |[Protein Generation Task] Protein novelty results after BLAST |