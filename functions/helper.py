import os
from pysanger import alignment

def align_row(row, input_folder, template_folder):
    """
    Aligns a sequencing file to a template file on a specified strand.

    Args:
        row (pandas.Series): A pandas Series object representing a row of data.
        input_folder (str): The path to the folder containing the .abi sequencing file.
        template_folder (str): The path to the folder containing the .dna, .xdna or .gb vector map templates.

    Returns:
        tuple: A tuple containing the sample name, template name and the alignment result.

    Raises:
        Exception: If an error occurs during the alignment process.

    """
    try:
        sample_name = row['Seq_File']
        template_name = row['Template_File']
        strand = row['strand']
        print(f"Aligning sequencing file `{sample_name}` to template `{template_name}` on strand `{strand}`.")
        alignment_result = alignment(
            abidata=os.path.join(input_folder, sample_name),
            template=os.path.join(template_folder, template_name),
            strand=strand
        )
        return sample_name,template_name, alignment_result
    
    except Exception as e:
        print(f"Error processing {row}: {e}")
        return row['Seq_File'], None
