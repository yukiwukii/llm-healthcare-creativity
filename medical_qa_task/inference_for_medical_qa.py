import ollama
import pandas as pd
import os
from tqdm import tqdm

# Define file paths
num = 20
input_file = f"mashqa/mashqa_part_{num}.csv"
final_output_file = f"mashqa/final_part_{num}.csv"
checkpoint_file = f"mashqa/checkpoint_part_{num}.csv"

# Ensure the 'mashqa' directory exists
if not os.path.exists('mashqa'):
    os.makedirs('mashqa')

# List of models to use
models = ['llama3.1:8b', 'gemma2:9b', 'mistral']

# Check for existing checkpoints
if os.path.exists(checkpoint_file):
    print(f"Resuming from checkpoint: {checkpoint_file}")
    df_checkpoint = pd.read_csv(checkpoint_file)
    df = pd.read_csv(input_file)

    # Merge the checkpoint data into the main dataframe
    df = df.merge(df_checkpoint, on="question", how="left", suffixes=("", "_checkpoint"))

    # Fill missing values in main columns with checkpoint data
    for model in models:
        column_name = f"{model.split(':')[0]}_answer"
        if column_name in df_checkpoint.columns:
            df[column_name] = df[column_name].combine_first(df_checkpoint[column_name])

    # Determine the start index based on completed rows in the checkpoint
    start_index = len(df_checkpoint.dropna(subset=[f"{model.split(':')[0]}_answer" for model in models]))
else:
    print("No checkpoint found. Starting from the beginning.")
    df = pd.read_csv(input_file)
    start_index = 0

# Initialize the prompt
prompt = '''You are a smart doctor. You are able to answer difficult medical questions.
    You are given a medical question. You must be straightforward in your answer.
    Do not give any links. Generate a complete answer in words. Answer with multiple answers.
    Each answer must be different, diverse, and numbered.
    If you generate same answers, or unfinished answers, you will be penalized.
    Do not use bold font.

    Example:
    Question:
    What questions might my doctor ask about my erectile dysfunction?
    Answer:
    1. The questions may include: Do you ever get an erection \n2. If you do, is it firm enough to have sex \n3. If you do start to have sex, do you then lose the erection \n4. Does it ever come back \n5. Can you get an erection by masturbation \n6. Do you ever wake up with an erection \n7. The doctor will ask if you smoke, how much alcohol you drink, and whether or not you use recreational drugs.'''

# Process the remaining rows
for index, row in tqdm(df.iloc[start_index:].iterrows(), total=len(df) - start_index, desc="Processing questions"):
    question = row["question"]
    for model in models:
        try:
            # Call the Ollama chat model
            response = ollama.chat(
                model=model,
                messages=[
                    {"role": "system", "content": prompt},
                    {"role": "user", "content": question}
                ]
            )
            # Update the corresponding column in the dataframe
            column_name = f"{model.split(':')[0]}_answer"
            df.loc[index, column_name] = response['message']['content']
        except Exception as e:
            # Handle exceptions and log them
            column_name = f"{model.split(':')[0]}_answer"
            df.loc[index, column_name] = f"Error: {str(e)}"
            with open("mashqa/error_log.txt", "a") as log_file:
                log_file.write(f"Row {index}, Model {model}: {str(e)}\n")

    # Save checkpoint every 10 rows
    if (index + 1) % 100 == 0:
        df.iloc[:index + 1].to_csv(checkpoint_file, index=False)

# Save the final dataset
df.to_csv(final_output_file, index=False)
print(f"Final dataset saved to {final_output_file}")
