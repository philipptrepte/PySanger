from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from snapgene_reader import snapgene_file_to_dict
from plotly.subplots import make_subplots
from functions.helper import abi_to_dict
from functions.helper import colorbar
from functions.helper import colorbar_plotly
from functions.helper import generate_consensusseq
from functions.helper import adjust_quality_scores

def alignment(abidata=None, template=None, strand=1):  
    """
    Perform sequence alignment between sequencing data and a reference sequence or a user-provided sequence.

    Parameters:
    - abidata (str): The path to the ABI file.
    - template (str or file path): A reference sequence in GenBank (.gb), DNA (.dna), or XDna (.xdna) format, or a user-provided DNA sequence (str).
    - strand (int): The strand to align the sequencing data to. Default is 1 (forward strand), -1 can be used for reverse strand alignment.

    Returns:
    - _abidata (dict): The extracted sequencing data from the ABI file from the abi_to_dict(filename) function.
    - asubject (str): The aligned subject sequence.
    - subject (str): The original subject sequence.
    - atemplate (str): The aligned template sequence.
    - sequence (str): The template sequence.
    - matches (list): A list indicating the match status of each base pair in the alignment. 1 for match, -1 for mismatch, 0 for not aligned.
    - ss (int): The number of gaps at the start of the subject sequence.
    - se (int): The number of gaps at the end of the subject sequence.
    - ts (int): The number of gaps at the start of the template sequence.
    - te (int): The number of gaps at the end of the template sequence.
    - avalues (list): The A channel values from the ABI file.
    - tvalues (list): The T channel values from the ABI file.
    - gvalues (list): The G channel values from the ABI file.
    - cvalues (list): The C channel values from the ABI file.
    - name (str): The name of the template sequence.
    - seq_filename (str): The filename of the ABI file.

    Raises:
    - ValueError: If the ABI file is invalid or the template is invalid.

    """
    if abidata is not None:
        _abidata, seq_filename = abi_to_dict(abidata)
    else:
        raise ValueError("Invalid ABI file. Please provide a valid ABI file.")
    
    avalues = _abidata["_channel"]["A"]
    tvalues = _abidata["_channel"]["T"]
    gvalues = _abidata["_channel"]["G"]
    cvalues = _abidata["_channel"]["C"]
    consensus_seq_set = generate_consensusseq(_abidata)
    
    if strand == 1: 
        subject = consensus_seq_set[0]
    
    elif strand == -1:
        subject = consensus_seq_set[1] 
        #avalues, tvalues, gvalues, cvalues = [list(reversed(x)) for x in (avalues, tvalues, gvalues, cvalues)]
        avalues, tvalues, gvalues, cvalues = (
            list(reversed(avalues)),
            list(reversed(tvalues)),
            list(reversed(gvalues)),
            list(reversed(cvalues))
        )
        _abidata["conf"] = list(reversed(_abidata["conf"]))

    if isinstance(template, str) and template.endswith((".gb", ".dna", ".xdna")):
        print(f"Reading vector map {template} as template.")
        if template.endswith(".gb"):
            """Reads a .gb file and transforms all sequences to upper case."""
            with open(template, "r") as file:
                sequence = SeqIO.read(file, "genbank")
                sequence = str(sequence.seq).upper()
                name = str(template)
                #return sequence.seq, name
                
        elif template.endswith(".dna"):
            """Reads a .dna file and transforms all sequences to upper case."""
            data = snapgene_file_to_dict(filepath=template)
            sequence = Seq(data["seq"])
            sequence = str(sequence).upper()
            name = str(template)
            #return sequence.seq, name
            
        elif template.endswith(".xdna"):
            """Reads a .xdna file and transforms all sequences to upper case."""
            with open(template, "rb") as file:
                sequence = next(XdnaIterator(file))
                sequence = str(sequence.seq).upper()
                name = str(template)
                #return sequence.seq, name    
                
        else:
            raise ValueError("Invalid file format. Please provide a .gb, .dna, or .xdna file.")
    
    elif type(template) == str:
        print("Reading user-provided sequence as template.")
        sequence = template.upper()
        name = "User-provided sequence"

    else:
        raise ValueError("Invalid template. Please provide a DNA string ('ATGC') or a file path.") 

    if sequence is not None:
        # Create a pairwise sequence aligner
        aligner = PairwiseAligner()
        
        # Set alignment parameters
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = 0
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -1
        aligner.end_open_gap_score = 0
        aligner.end_extend_gap_score = 0

        if name == "User-provided sequence":
            print(f"Aligning sequencing file `{seq_filename}` ({len(subject)} bp) to `{name}` ({len(sequence)} bp).")
        
        elif name != "User-provided sequence":
            print(f"Aligning sequencing file `{seq_filename}` ({len(subject)} bp) to provided reference sequence `{name}` ({len(sequence)} bp).")
        
        alignments  = aligner.align(subject, sequence)
        #print(alignments[0])
        print(f"Alignment score: {alignments[0].score}")
        
        atemplate   = alignments[0][1].upper()
        asubject    = alignments[0][0].upper()
        

        ts = 0 
        ss = 0 
        for t, s in zip(atemplate, asubject):
            if t != "-" and s == "-": 
                ts += 1
            
            elif t == "-" and s != "-": 
                ss += 1  
            
            elif t == "-" and s == "-":
                ts += 1
                ss += 1

            else:
                break 

        te = 0 
        se = 0 
        for t, s in zip(atemplate[::-1], asubject[::-1]):
            if t != "-" and s == "-": 
                te += 1
            elif t == "-" and s != "-": 
                se += 1  
            elif t == "-" and s == "-":
                te += 1
                se += 1
            else:
                break 
       
        pos = 0 
        tpos = 0
        matches = []
        new_avalues = []
        new_tvalues = [] 
        new_gvalues = []
        new_cvalues = [] 
        
        for t, s in zip(atemplate, asubject):
            if s == "-":
                new_avalues.extend([0,0,0,0,0]) 
                new_tvalues.extend([0,0,0,0,0]) 
                new_gvalues.extend([0,0,0,0,0]) 
                new_cvalues.extend([0,0,0,0,0])     
            else:
                new_avalues.append(avalues[5*pos]) 
                new_avalues.append(avalues[5*pos+1]) 
                new_avalues.append(avalues[5*pos+2])
                new_avalues.append(avalues[5*pos+3]) 
                new_avalues.append(avalues[5*pos+4])

                new_tvalues.append(tvalues[5*pos])
                new_tvalues.append(tvalues[5*pos+1]) 
                new_tvalues.append(tvalues[5*pos+2])
                new_tvalues.append(tvalues[5*pos+3]) 
                new_tvalues.append(tvalues[5*pos+4])
                
                new_gvalues.append(gvalues[5*pos])
                new_gvalues.append(gvalues[5*pos+1]) 
                new_gvalues.append(gvalues[5*pos+2])
                new_gvalues.append(gvalues[5*pos+3]) 
                new_gvalues.append(gvalues[5*pos+4])

                new_cvalues.append(cvalues[5*pos])
                new_cvalues.append(cvalues[5*pos+1])
                new_cvalues.append(cvalues[5*pos+2])
                new_cvalues.append(cvalues[5*pos+3]) 
                new_cvalues.append(cvalues[5*pos+4])

                pos += 1
            tpos += 1

            if t == s: 
                matches.append(1) 
            elif len(subject) > len(sequence) and s != t and ss < tpos < ss+len(sequence):
                matches.append(-1)
            elif len(sequence) > len(subject) and s != t and ts < tpos < ts+len(subject):
                matches.append(-1)
            else:
                matches.append(0) 
            
            match_count = matches.count(1)
            mismatch_count = matches.count(-1)
            not_aligned_count = matches.count(0)
        print(f"{match_count} matches found between base pairs {ts} and {ts+len(subject)} in '{name}'. {mismatch_count} mismatches found. {not_aligned_count} base pairs not aligned.")
        
        avalues = new_avalues 
        tvalues = new_tvalues
        gvalues = new_gvalues
        cvalues = new_cvalues
    return _abidata, asubject, subject, atemplate, sequence, matches, ss, se, ts, te, avalues, tvalues, gvalues, cvalues, name, seq_filename, strand

def visualize(alignment, region="aligned", fontsize=2, fig_width=None):
    """
    Visualizes the alignment data.

    Args:
        alignment (tuple): A tuple containing the alignment data from the alignment(abidata, template, strand) function.
        region (str, optional): The region of the alignment to visualize. Can be either "all" or "aligned". Defaults to "aligned".
        fontsize (float or int, optional): The font size for the visualization. Defaults to 2.
        fig_width (float or int, optional): The width of the figure. Defaults to None. If None, the width is automatically determined based on the sequence length.

    Returns:
        matplotlib.figure.Figure: The generated figure.
    """
    _abidata, asubject, subject, atemplate, sequence, matches, ss, se, ts, te, avalues, tvalues, gvalues, cvalues, name, seq_filename, strand = alignment
    plt.close("all")

    assert isinstance(fontsize, (float, int)), f"fontsize must be a float or integer, not '{fontsize}'."
    assert isinstance(fig_width, (float, int)) or fig_width is None, f"Figure width must be a float or integer, not '{fig_width}'."

    if region == "all":
        if fig_width is None:
            fig_width = 50
        tasubject  = asubject
        tatemplate = atemplate
        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(fig_width, 3), gridspec_kw={'height_ratios': [0.5, 3, 0.5]})

    elif region == "aligned":
        if fig_width is None:
            raise ValueError("Figure width must be provided when visualizing the aligned region.")
        if ts == 0 and te == 0:
            tasubject  = asubject[ss:len(asubject) - se]
            tatemplate = atemplate[ss:len(atemplate) - se]
            matches   = matches[ss:len(asubject) - se]
        else:
            tasubject  = asubject[ts:len(asubject) - te]
            tatemplate = atemplate[ts:len(atemplate) - te]
            matches   = matches[ts:len(asubject) - te]
        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(fig_width, 3), gridspec_kw={'height_ratios': [0.5, 3, 0.5]})

    else:
        raise ValueError("Invalid region. Please provide 'all' or 'aligned'.")
    
    ax = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]
    
    # Determine tick spacing based on the length of the sequences
    if region == "all":
        tlength = len(sequence)
    elif region == "aligned":
        tlength = len(tatemplate)
    if tlength < 100:
        ttick_space = 10
    elif tlength < 200:
        ttick_space = 20
    elif tlength < 500:
        ttick_space = 50
    elif tlength < 1000:
        ttick_space = 50
    elif tlength < 2500:
        ttick_space = 100    
    else:
        ttick_space = 250
    if ttick_space > len(sequence):
        if len(sequence) < 100:
            ttick_space = 10
        elif len(sequence) < 200:
            ttick_space = 20
        elif len(sequence) < 500:
            ttick_space = 50
        elif len(sequence) < 1000:
            ttick_space = 50
        elif len(sequence) < 2500:
            ttick_space = 100    

    if tatemplate is not None:
        ax = colorbar(ax, tatemplate, matches=matches, char=True, fontsize=fontsize)
        tnum = 0
        if region == "all":
            if ts != 0:
                ss = ts
                ts = 0
            elif ts == 0:
                ts = ss
                ss = 0
            tticks      = [] 
            tticklabels = [] 
            for letter in tatemplate:
                if (tnum) % ttick_space == 0: 
                    tticks.append(ts + tnum) 
                    tticklabels.append(str(tnum))
                    tnum += 1
                if letter == "-":
                    pass 
                else:
                    tnum += 1 
        elif region == "aligned":
            tticks      = [0.5]
            tticklabels = [ts]
            for letter in tatemplate:
                if (tnum+1) % ttick_space == 0: 
                    tticks.append(tnum+0.5) 
                    tticklabels.append(str(ts+tnum+1))
                    tnum += 1
                if letter == "-":
                    pass 
                else:
                    tnum += 1 
        ax.xaxis.tick_top() 
        ax.set_xticks(tticks) 
        ax.set_xticklabels(tticklabels) 
        
        ax.set_ylabel(name, rotation=0, va="center", ha="right") 
    else:
        raise ValueError("Invalid template. No template sequence provided. Please check alignment.")
    
    if region == "all":
        quality = adjust_quality_scores(tasubject, quality_scores = _abidata["conf"])
    elif region == "aligned":
        quality = adjust_quality_scores(tasubject, quality_scores = _abidata["conf"])
    
    if region == "all":
        sns.barplot(x=list(range(len(tasubject))), y=quality, color="#F7F7FF", edgecolor="#BBBBBB", linewidth=0.5, ax=ax2, width = 1.0)
        ax2.set_xlim(-0.5, len(tatemplate) - 0.5) 
    elif region == "aligned":
        if ts == 0 and te == 0:
            sns.barplot(x=list(range(len(tasubject))), y=quality[ss:-se], color="#F7F7FF", edgecolor="#BBBBBB", linewidth=0.5, ax=ax2, width = 1.0)
        elif ts != 0 and te != 0:
            sns.barplot(x=list(range(len(tasubject))), y=quality, color="#F7F7FF", edgecolor="#BBBBBB", linewidth=0.5, ax=ax2, width = 1.0)
        ax2.set_xlim(-0.5, len(tasubject) - 0.5) 
    
    ax2X2 = ax2.twinx()

    if region == "all":
        min_length = min(len(tvalues), len(avalues), len(gvalues), len(cvalues))
        positions = np.linspace(-0.5, (min_length / 5) - 0.5, min_length).tolist()
        sns.lineplot(x=positions, y=tvalues, color="#FC58FE", lw=0.5, ax=ax2X2)
        sns.lineplot(x=positions, y=avalues, color="#33CC33", lw=0.5, ax=ax2X2)
        sns.lineplot(x=positions, y=gvalues, color="#303030", lw=0.5, ax=ax2X2)
        sns.lineplot(x=positions, y=cvalues, color="#395CC5", lw=0.5, ax=ax2X2)
        combined_values = tvalues + avalues + gvalues + cvalues

    elif region == "aligned":
        if ts == 0 and te == 0:
            _te = -5 * se
            ts = ss
        else:
            _te = -5 * te
        # Ensure the lengths of the arrays are consistent
        min_length = min(len(tvalues[5 * ts:_te]), len(avalues[5 * ts:_te]), len(gvalues[5 * ts:_te]), len(cvalues[5 * ts:_te]))
        positions = np.linspace(-0.5, (min_length / 5) - 0.5, min_length).tolist()
        sns.lineplot(x=positions, y=tvalues[5 * ts:5 * ts + min_length], color="#FC58FE", lw=0.5, ax=ax2X2)
        sns.lineplot(x=positions, y=avalues[5 * ts:5 * ts + min_length], color="#33CC33", lw=0.5, ax=ax2X2)
        sns.lineplot(x=positions, y=gvalues[5 * ts:5 * ts + min_length], color="#303030", lw=0.5, ax=ax2X2)
        sns.lineplot(x=positions, y=cvalues[5 * ts:5 * ts + min_length], color="#395CC5", lw=0.5, ax=ax2X2)

        #combined_values = tvalues[5 * ss:5 * ss + min_length] + avalues[5 * ss:5 * ss + min_length] + gvalues[5 * ss:5 * ss + min_length] + cvalues[5 * ss:5 * ss + min_length]
        combined_values = tvalues[5 * ts:5 * ts + min_length] + avalues[5 * ts:5 * ts + min_length] + gvalues[5 * ts:5 * ts + min_length] + cvalues[5 * ts:5 * ts + min_length]

    ax2.set_xticklabels("")
    ax2.set_xticks([])
    ax2.set_ylim(0, 1.01 * max(quality))
    ax2.set_ylabel("Quality")
    ax2X2.set_ylim(min(combined_values), 1.01 * max(combined_values))
    ax2X2.set_ylabel("Peak height")

    ax2.spines["top"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2X2.spines["top"].set_visible(False)
    ax2X2.spines["bottom"].set_visible(False)
    ax2X2.spines["left"].set_visible(False)

    ax3 = colorbar(ax3, tasubject, matches=matches, char=True, fontsize=fontsize)

    if region == "all":
        slength = len(subject)
    elif region == "aligned":
        slength = len(tasubject)
    if slength < 100:
        stick_space = 10
    elif slength < 200:
        stick_space = 20
    elif slength < 500:
        stick_space = 50
    elif slength < 1000:
        stick_space = 100
    elif slength < 2500:
        stick_space = 250    
    else:
        stick_space = 500
    snum = 0 
    if region == "all":
        sticks      = [ss] 
        sticklabels = [ss] 
        for letter in tasubject:
            if (ss+snum+1) % stick_space == 0: 
                sticks.append(ss+snum+1) 
                sticklabels.append(str(ss+snum+1))
                snum += 1
            if letter == "-":
                pass 
            else:
                snum += 1 
    elif region == 'aligned':
        if strand == 1:
            if ts != 0 and te != 0:
                ts = ss
            sticks      = [0.5] 
            sticklabels = [ss]
            for letter in tasubject:
                if (ts + snum+1) % stick_space == 0: 
                    sticks.append(snum+0.5) 
                    sticklabels.append(str(ss+snum+1))
                if letter == "-":
                    pass 
                else:
                    snum += 1
        if strand == -1:
            snum = len(subject)
            ts = len(subject)
            sticks = [snum]
            sticklabels = [ss]
            for letter in tasubject[::-1]:
                if (snum) % stick_space == 0: 
                    sticks.append(snum-0.5) 
                    sticklabels.append(str(ts-snum+1))
                if letter == "-":
                    pass 
                else:
                    snum -= 1
    

    ax3.set_xticks(sticks) 
    ax3.set_xticklabels(sticklabels) 
    ax3.set_ylabel(seq_filename, rotation=0, va="top", ha="right") 

    # Adjust the spacing between subplots
    plt.subplots_adjust(hspace=0, wspace=0.3)  # Adjust hspace and wspace as needed

    #plt.tight_layout() # Adjust layout to ensure subplots fit into the figure area
    return fig

def visualize_plotly(alignment, region="aligned", fontsize=10):
    """"
    Generates a Plotly figure to visualize sequence alignment data.

    Parameters:
    - alignment (tuple): A tuple containing the alignment data from the alignment(abidata, template, strand) function.
    - region (str): The region of the alignment to visualize. Can be either "all" or "aligned". Defaults to "aligned".
    - fontsize (int): The font size for the plot annotations. Defaults to 10.

    Returns:
    - fig (plotly.graph_objects.Figure): The generated Plotly figure.

    Raises:
    - ValueError: If an invalid region is provided. Valid options are "all" or "aligned".
    - ValueError: If no template sequence is provided in the alignment data.
    """

    _abidata, asubject, subject, atemplate, sequence, matches, ss, se, ts, te, avalues, tvalues, gvalues, cvalues, name, seq_filename, strand = alignment

    if region == "all":
        tasubject  = asubject
        tatemplate = atemplate
        fig_height = 600
    elif region == "aligned":
        if(min(len(sequence), len(subject)) < 100):
            fig_width = 500
        else:
            fig_width = 2000
        if ts == 0 and te == 0:
            tasubject  = asubject[ss:len(asubject) - se]
            tatemplate = atemplate[ss:len(atemplate) - se]
            matches   = matches[ss:len(asubject) - se]
        elif ts != 0 and te != 0:
            tasubject  = asubject[ts:len(asubject) - te]
            tatemplate = atemplate[ts:len(atemplate) - te]
            matches   = matches[ts:len(asubject) - te]
        fig_height = 600
    else:
        raise ValueError("Invalid region. Please provide 'all' or 'aligned'.")

    fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.02, row_heights=[0.2, 0.6, 0.2], specs=[[{"secondary_y": False}], [{"secondary_y": True}], [{"secondary_y": False}]])

    # Determine tick spacing based on the length of the sequences
    if region == "all":
        tlength = len(sequence)
    elif region == "aligned":
        tlength = len(tatemplate)
    if tlength < 100:
        ttick_space = 10
    elif tlength < 200:
        ttick_space = 20
    elif tlength < 500:
        ttick_space = 50
    elif tlength < 1000:
        ttick_space = 50
    elif tlength < 2500:
        ttick_space = 100    
    else:
        ttick_space = 250
    if ttick_space > len(sequence):
        if len(sequence) < 100:
            ttick_space = 10
        elif len(sequence) < 200:
            ttick_space = 20
        elif len(sequence) < 500:
            ttick_space = 50
        elif len(sequence) < 1000:
            ttick_space = 50
        elif len(sequence) < 2500:
            ttick_space = 100    

    if tatemplate is not None:
        ax, annotations1 = colorbar_plotly(tatemplate, matches=matches, char=True, fontsize=fontsize)
        for annotation in annotations1:
            annotation.update(dict(xref='x1', yref='y1'))  # Ensure annotations reference the first subplot
            fig.add_annotation(annotation)
        tnum = 0
        if region == "all":
            if ts != 0:
                ss = ts
                ts = 0
            elif ts == 0:
                ts = ss
                ss = 0
            tticks      = [] 
            tticklabels = [] 
            for letter in tatemplate:
                if (tnum) % ttick_space == 0: 
                    tticks.append(ts + tnum) 
                    tticklabels.append(str(tnum))
                    tnum += 1
                if letter == "-":
                    pass 
                else:
                    tnum += 1 
        elif region == "aligned":
            tticks      = [0.5]
            tticklabels = [ts]
            for letter in tatemplate:
                if (tnum+1) % ttick_space == 0: 
                    tticks.append(tnum+0.5) 
                    tticklabels.append(str(ts+tnum+1))
                    tnum += 1
                if letter == "-":
                    pass 
                else:
                    tnum += 1 
        
        for trace in ax.data:
            fig.add_trace(trace, row=1, col=1)
        
        fig.update_layout(
            xaxis1=dict(showgrid=False, zeroline=False, tickvals=tticks, ticktext=tticklabels, showticklabels = True, side='top'),
            yaxis1=dict(showticklabels=False, showgrid=False, zeroline=False, title_text=name),
            margin=dict(l=0, r=0, t=0, b=0),
            plot_bgcolor='rgba(0,0,0,0)',
            paper_bgcolor='rgba(0,0,0,0)'
        )
                
    else:
        raise ValueError("Invalid template. No template sequence provided. Please check alignment.")
    
    if region == "all":
        quality = adjust_quality_scores(tasubject, quality_scores = _abidata["conf"])
    elif region == "aligned":
        quality = adjust_quality_scores(tasubject, quality_scores = _abidata["conf"])
    
    if region == "all":
        fig.add_trace(go.Bar(x=list(range(len(tasubject))), y=quality, marker_color="#F7F7FF", marker_line_color="#BBBBBB", marker_line_width=0.5), row=2, col=1)
    elif region == "aligned":
        if ts == 0 and te == 0:
            fig.add_trace(go.Bar(x=list(range(len(tasubject))), y=quality[ss:-se], marker_color="#F7F7FF", marker_line_color="#BBBBBB", marker_line_width=0.5), row=2, col=1)
        elif ts != 0 and te != 0:
            fig.add_trace(go.Bar(x=list(range(len(tasubject))), y=quality, marker_color="#F7F7FF", marker_line_color="#BBBBBB", marker_line_width=0.5), row=2, col=1)
    
    if region == "all":
        min_length = min(len(tvalues), len(avalues), len(gvalues), len(cvalues))
        positions = np.linspace(-0.5, (min_length / 5) - 0.5, min_length).tolist()
        fig.add_trace(go.Scatter(x=positions, y=tvalues, mode='lines', line=dict(color="#FC58FE", width=1)), row=2, col=1, secondary_y=True)
        fig.add_trace(go.Scatter(x=positions, y=avalues, mode='lines', line=dict(color="#33CC33", width=1)), row=2, col=1, secondary_y=True)
        fig.add_trace(go.Scatter(x=positions, y=gvalues, mode='lines', line=dict(color="#303030", width=1)), row=2, col=1, secondary_y=True)
        fig.add_trace(go.Scatter(x=positions, y=cvalues, mode='lines', line=dict(color="#395CC5", width=1)), row=2, col=1, secondary_y=True)
        combined_values = tvalues + avalues + gvalues + cvalues
    elif region == "aligned":
        if ts == 0 and te == 0:
            _te = -5 * se
            ts = ss
        elif ts !=0 and te != 0:
            _te = -5 * te
        min_length = min(len(tvalues[5 * ts:_te]), len(avalues[5 * ts:_te]), len(gvalues[5 * ts:_te]), len(cvalues[5 * ts:_te]))
        positions = np.linspace(-0.5, (min_length / 5) - 0.5, min_length).tolist()
        fig.add_trace(go.Scatter(x=positions, y=tvalues[5 * ts:5 * ts + min_length], mode='lines', line=dict(color="#FC58FE", width=1)), row=2, col=1, secondary_y=True)
        fig.add_trace(go.Scatter(x=positions, y=avalues[5 * ts:5 * ts + min_length], mode='lines', line=dict(color="#33CC33", width=1)), row=2, col=1, secondary_y=True)
        fig.add_trace(go.Scatter(x=positions, y=gvalues[5 * ts:5 * ts + min_length], mode='lines', line=dict(color="#303030", width=1)), row=2, col=1, secondary_y=True)
        fig.add_trace(go.Scatter(x=positions, y=cvalues[5 * ts:5 * ts + min_length], mode='lines', line=dict(color="#395CC5", width=1)), row=2, col=1, secondary_y=True)
        combined_values = tvalues[5 * ss:5 * ss + min_length] + avalues[5 * ss:5 * ss + min_length] + gvalues[5 * ss:5 * ss + min_length] + cvalues[5 * ss:5 * ss + min_length]

    fig.update_yaxes(title_text="Quality", row=2, col=1)
    fig.update_yaxes(title_text="Peak height", row=2, col=1, secondary_y=True)

    if region == "all":
        slength = len(subject)
    elif region == "aligned":
        slength = len(tasubject)
    if slength < 100:
        stick_space = 10
    elif slength < 200:
        stick_space = 20
    elif slength < 500:
        stick_space = 50
    elif slength < 1000:
        stick_space = 100
    elif slength < 2500:
        stick_space = 250    
    else:
        stick_space = 500
    ax3, annotations3 = colorbar_plotly(tasubject, matches=matches, char=True, fontsize=fontsize)
    for annotation in annotations3:
        annotation.update(dict(xref='x3', yref='y4'))  # Update xref and yref to refer to the third subplot
        fig.add_annotation(annotation)
        
    snum = 0 
    if region == "all":
        sticks      = [ss] 
        sticklabels = [ss] 
        for letter in tasubject:
            if (ss+snum+1) % stick_space == 0: 
                sticks.append(ss+snum+1) 
                sticklabels.append(str(ss+snum+1))
                snum += 1
            if letter == "-":
                pass 
            else:
                snum += 1 
    elif region == 'aligned':
        if ts != 0 and te != 0:
            ts = ss
        sticks      = [0.5] 
        sticklabels = [ss]
        for letter in tasubject:
            if (ts + snum+1) % stick_space == 0: 
                sticks.append(snum+0.5) 
                sticklabels.append(str(ss+snum+1))
            if letter == "-":
                pass 
            else:
                snum += 1

    for trace in ax3.data:
            fig.add_trace(trace, row=3, col=1)
    fig.update_layout(
        xaxis3=dict(showgrid=False, zeroline=False, tickvals=sticks, ticktext=sticklabels, showticklabels = True, side='bottom'),
        yaxis4=dict(showticklabels=False, showgrid=False, zeroline=False, title_text=seq_filename),
        margin=dict(l=0, r=0, t=0, b=0),
        plot_bgcolor='rgba(0,0,0,0)',
        paper_bgcolor='rgba(0,0,0,0)'
    )    

    fig.update_layout(height=fig_height, showlegend=False)
    return fig

if __name__ == "__main__":

    align        = alignment(abidata="examples/BE MAFB5.ab1", template="AGCCGGCTGGCTGCAGGCGT")
    fig          = visualize(align, region = "all", fontsize = 5)
    fig.show()
    fig2         = visualize(align, region = "aligned", fontsize = 10)
    fig2.show()

    align2        = alignment(abidata="seq_results/QPSQ0664-CMV-for.ab1", template="templates/QPPL0052_pcDNA3.1_mCitrine-C1-GW.dna", strand = 1)
    fig2          = visualize(align2, region="all")
    fig2.show()
    
    align3        = alignment(abidata="seq_results/QPSQ0665-BGHR.ab1", template="templates/QPPL0052_pcDNA3.1_mCitrine-C1-GW.dna", strand = -1)
    fig3          = visualize(align3, region="aligned")
    fig3.show()

    fig4 = visualize_plotly(align2, region="aligned")
    fig4.show()