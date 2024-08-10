import os
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import plotly.graph_objects as go

def abi_to_dict(filename):
    """
        Convert an ABI file to a dictionary containing sequencing data.

        This function reads an ABI file and extracts the sequencing data, including confidence values and channel intensities,
        and stores it in a dictionary format. The dictionary structure includes the following keys:
        
        - 'conf': A list of confidence values for each base position.
        - 'channel': A nested dictionary containing intensity values for each base (A, T, G, C) at each position.
        - '_channel': A nested dictionary containing intensity values for each base (A, T, G, C) at each position, including
                                    neighboring positions.

        Parameters:
        - filename (str): The path to the ABI file.

        Returns:
        - abi_data (dict): The dictionary containing the extracted sequencing data.
        - filename (str): The input filename.

        """
    
    record   = SeqIO.read(filename,'abi')
    abi_data = {"conf":[],
                "channel":{"A":[],
                           "T":[],
                           "G":[],
                           "C":[],
                          },
                "_channel":{"A":[],
                            "T":[],
                            "G":[],
                            "C":[],
                          }
                }
    for i, (pos, conf) in enumerate(zip(record.annotations['abif_raw']["PLOC1"], record.annotations['abif_raw']["PCON1"])):
        if pos > 4 and pos < len(record.annotations['abif_raw']["DATA9"])-5: 
            abi_data["conf"].append(conf)
            abi_data["channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])

            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-5])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+5])
            
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-5])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+5])
            
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-5])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+5])

            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-5])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+5])

    return abi_data, filename

def generate_consensusseq(abidata):
    """
    Generate a consensus sequence based on the given abidata.

    Args:
        abidata (dict): A dictionary containing the abidata.

    Returns:
        tuple: A tuple containing two strings. The first string is the consensus sequence,
               and the second string is the reverse complement of the consensus sequence.
    """
    consensus_seq = "" 
    _atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}
    
    for values in zip(abidata["channel"]["A"], abidata["channel"]["T"], abidata["channel"]["G"], abidata["channel"]["C"]):
        consensus_seq += _atgc_dict[values.index(max(values))]
     
    #return (consensus_seq, consensus_seq.translate(str.maketrans("ATGC","TACG"))[::-1])
    return (Seq(consensus_seq), Seq(consensus_seq).reverse_complement())

def generate_pwm(abidata):
    """
    Generate a position weight matrix (PWM) based on the given abidata.

    Args:
        abidata (tuple): A tuple containing the abidata and seq_filename.

    Returns:
        pandas.DataFrame: The generated PWM as a pandas DataFrame.

    """
    _abidata, seq_filename = abidata
    pwm = {"A":[], "T":[], "G":[], "C":[]} 
    for values in zip(_abidata["channel"]["A"], _abidata["channel"]["T"], _abidata["channel"]["G"], _abidata["channel"]["C"]):
        v = 100000 / (sum(values)+1) 
        new_values = (v*values[0], v*values[1], v*values[2], v*values[3])
        new_values = list(map(int, new_values))
        
        while sum(new_values) < 100000:
            for i in range(len(new_values)):
                new_values[i] += 1
                if sum(new_values) == 100000:
                    break 
        
        pwm["A"].append(new_values[0])
        pwm["T"].append(new_values[1])
        pwm["G"].append(new_values[2])
        pwm["C"].append(new_values[3])
    
    pwm=pd.DataFrame(pwm)
    return pwm

def colorbar(ax, ref, matches=None, char=True, fontsize=2):
    """
    Add a color bar to the given axes object.

    Parameters:
    - ax (matplotlib.axes._axes.Axes): The axes object to which the color bar will be added.
    - ref (list): A list of reference values.
    - matches (list, optional): A list of match values. If provided, the color bar will be customized based on the matches.
    - char (bool, optional): Whether to display characters on the color bar. Default is True.
    - fontsize (int, optional): The font size of the characters on the color bar. Default is 2.

    Returns:
    - matplotlib.axes.Axes: The modified axes object.

    """
    bars = ax.bar(list(range(len(ref))), [0.9] * (len(ref)), width=1.0, edgecolor="#BBBBBB", linewidth=0.5, align="edge",bottom=0.05)
    ax.set_xlim(0,len(ref))
    ax.set_ylim(0,1.00)
    p = 0
    if matches is None:
        for bar, c in zip(bars,ref):
            bar.set_facecolor("w")
            if char == True:
                ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
            p += 1
    
    else:
        for m, bar, c in zip(matches, bars, ref):
            if m == -1:
                bar.set_alpha(0.5)
                bar.set_facecolor("r")
                bar.set_edgecolor("#BBBBBB")
                bar.set_linewidth(0.5)
            else:
                bar.set_facecolor("w")
                bar.set_edgecolor("#BBBBBB")
                bar.set_linewidth(0.5)
            
            if char == True:
                ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
            p += 1
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    return ax

def colorbar_plotly(ref, matches=None, char=True, fontsize=12):
    """
    Generate a color bar plot using Plotly.

    Args:
        ref (str): The reference sequence.
        matches (list, optional): A list of matches. Default is None.
        char (bool, optional): Whether to display characters on the color bar. Default is True.
        fontsize (int, optional): The font size of the characters. Default is 12.

    Returns:
        tuple: A tuple containing the Plotly figure object and a list of annotations.
    """
    bars = []
    annotations = []
    p = 0

    if matches is None:
        for c in ref:
            bars.append(go.Bar(
                x=[p],
                y=[0.9],
                width=1.0,
                marker=dict(color='white', line=dict(color='#BBBBBB', width=0.5)),
                showlegend=False
            ))
            if char:
                annotations.append(dict(
                    x=p,
                    y=0.45,
                    text=c,
                    showarrow=False,
                    font=dict(size=fontsize),
                    xanchor='center',
                    yanchor='middle'
                ))
            p += 1
    else:
        for m, c in zip(matches, ref):
            color = 'red' if m == -1 else 'white'
            alpha = 0.5 if m == -1 else 1.0
            bars.append(go.Bar(
                x=[p],
                y=[0.9],
                width=1.0,
                marker=dict(color=color, opacity=alpha, line=dict(color='#BBBBBB', width=0.5)),
                showlegend=False
            ))
            if char:
                annotations.append(dict(
                    x=p,
                    y=0.45,
                    text=c,
                    showarrow=False,
                    font=dict(size=fontsize),
                    xanchor='center',
                    yanchor='middle'
                ))
            p += 1

    fig = go.Figure(data=bars)

    return fig, annotations

def adjust_quality_scores(tasubject, quality_scores):
    """
    Adjusts the quality scores by inserting placeholder values at the positions of insertions in the alignment.

    Args:
        tasubject (str): The alignment sequence where insertions are represented by '-'.
        quality_scores (list): The list of quality scores corresponding to each base in the alignment.

    Returns:
        list: The adjusted quality scores with placeholder values inserted at the positions of insertions.
    """
    # Identify the positions of the insertions in the alignment
    insertion_positions = [i for i, base in enumerate(tasubject) if base == '-']
    
    # Adjust the quality scores by inserting placeholder values at the insertion positions
    adjusted_quality_scores = list(quality_scores)
    for pos in insertion_positions:
        adjusted_quality_scores.insert(pos, -1)  # Insert None or a specific low-quality score
    
    return adjusted_quality_scores

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
    from pysanger import alignment
    
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