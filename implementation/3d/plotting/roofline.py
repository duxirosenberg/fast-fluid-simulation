import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
from labellines import labelLine, labelLines
import numpy as np
import pandas as pd





darkgreen = (0.1,0.3,0.1,0.99)
# violet = # (0.2,0,0.3,0.99)



rc('font',**{'family':'serif','serif':['Helvetica']})
rc('text', usetex=True)
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Helvetica']
plt.rcParams['figure.figsize'] = [16, 9]



xmax = 64
ymax = 64
xmin = 1/64
ymin = 1/64
x_ticks = [1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32]
x_ticksL = ["1/16", "1/8", "1/4", "1/2", "1", "2", "4", "8", "16", "32"]
y_ticks = [1/16, 1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32]
y_ticksL = ["1/16", "1/8", "1/4", "1/2", "1", "2", "4", "8", "16", "32"]

# avx2 double, int
floats_per_vec = 4
ints_per_vec = 8


def get_data(ds, bc, L1, L2, L3):
    csv_file_path = "TimingData_roofline.csv"
    df = pd.read_csv(csv_file_path)
    df = df.loc[(df["function"]==bc)]
    df = df.loc[(df['DIRECTION_SIZE'] == ds)]
    nx = df["NX"].values
    iops = df["iops"].values
    flops = df["flops"].values
    ts = df["TIMESTEPS"].values
    cycles = df["cycles"].values

    bread = df["bytes_read"].values
    bwrite = df["bytes_write"].values
    
    bytes = bread + bwrite
    fOI = flops/bytes
    iOI = iops/bytes 
    fOIr = flops/bread
    iOIr = iops/bread
    fOIw = flops/bwrite
    iOIw = iops/bwrite 
    fP = ts*flops/cycles
    iP = ts*iops/cycles

    colors = np.where(bytes<L3, "purple", "darkviolet")
    colors = np.where(bytes<L2, "indigo", colors)
    colors = np.where(bytes<L1, "navy", colors)

    return fOI,iOI,fOIr,iOIr,fOIw,iOIw,fP,iP,nx,colors



# def get_colors


    # Peak CPU performance (ops/cycle)
    # Peak CPU performance Vectorized (ops/cycle)
    # Performance: Flops/cycle = flops/byte * bytes/cycles corresponding to the given intensities
    # Operational intensities (flops/byte)

def roofline(  
                peak_performance_flops,
                peak_performance_flops_fma,
                peak_performance_iops,
                L1_bw,
                L2_bw,
                L3_bw,
                mem_bw,
                fq,
                L1_s, 
                L2_s,
                L3_s,
                dType,
                fnames,
                op_types,
                byte_types,
                savefile,
                mycols,
                zen=False
                ):
    ax  = plt.gca()
    fig = plt.gcf()




    plt.loglog()
    plt.axis([xmin,xmax,ymin,ymax])
    plt.xticks(x_ticks, x_ticksL)
    plt.yticks(y_ticks, y_ticksL)
    plt.minorticks_off()
    plt.grid(True)

    plt.xlabel('Operational Intensity \n [iops/byte], [flops/byte]')
    plt.ylabel('Performance \n [iops/cycle], [flops/cycle]', rotation=0)
    ax.yaxis.set_label_coords(0.1, 1.0)
    ax.xaxis.set_label_coords(0.8, -0.1)
    plt.title('Roofline Plot - '+str(dType)+"Directions") #, x=0.85, y=0.9)






    intensity = np.logspace(-5, 1, num=1000) 

    peak_performance_flops_v = peak_performance_flops*floats_per_vec
    peak_performance_flops_v = peak_performance_flops*floats_per_vec
    peak_performance_flops_fma_v = peak_performance_flops_fma*floats_per_vec
    peak_performance_iops_v = peak_performance_iops*ints_per_vec



    # intel  amd slightly built different
    bw_performance_l1_read =  intensity * L1_bw[0]
    bw_performance_l1_write =  intensity * L1_bw[1]
    bw_performance_l2_rw =  intensity * (L2_bw)
    bw_performance_l3 =  intensity * L3_bw[0]
    bw_performance_mem = intensity * mem_bw




    ax.axhline(y=peak_performance_flops,        color=darkgreen, linestyle='--', label='peak Performance \n flop scalar')
    ax.axhline(y=peak_performance_flops_v,      color=darkgreen, linestyle='--', label='peak Performance \n flop vectorized')
    ax.axhline(y=peak_performance_flops_fma,    color=darkgreen, linestyle='--', label='peak Performance \n flop w.fma scalar')
    ax.axhline(y=peak_performance_flops_fma_v,  color=darkgreen, linestyle='--', label='peak Performance \n flop w. fma vectorized')
    ax.axhline(y=peak_performance_iops,         color=darkgreen, linestyle='--', label='peak Performance \n iop scalar')
    ax.axhline(y=peak_performance_iops_v,       color=darkgreen, linestyle='--', label='peak Performance \n iop vectorized')

    i = 0

    lines = ax.get_lines()
    labelLine(lines[i], 30, align=True, fontsize=9, yoffset=0.0);    i+=1
    labelLine(lines[i], 30, align=True, fontsize=9, yoffset=0.0);    i+=1
    labelLine(lines[i], 10, align=True, fontsize=9, yoffset=0.0);    i+=1
    labelLine(lines[i], 30, align=True, fontsize=9, yoffset=0.0);    i+=1
    labelLine(lines[i], 30, align=True, fontsize=9, yoffset=0.0);    i+=1
    labelLine(lines[i], 30, align=True, fontsize=9, yoffset=0.0);    i+=1


    if(zen):
        ax.plot(intensity, bw_performance_l1_read,  label="bound based on \n L1,2,3 read bw, write bw" , color ="navy", linestyle=':')
        lines = ax.get_lines();  labelLine(lines[i], 1/40, align=True, fontsize=9, yoffset=0.0) ;    i+=1
        ax.plot(intensity, bw_performance_l2_rw,    label="bound based on \n L1,2,3 read+write bandwith", color ="navy", linestyle=':')
        lines = ax.get_lines();  labelLine(lines[i], 1/40, align=True, fontsize=9, yoffset=0.0)  ;   i+=1
        
    else:
        ax.plot(intensity, bw_performance_l1_read,  label="bound based on \n L1 read bandwith" , color ="navy", linestyle=':')
        lines = ax.get_lines();  labelLine(lines[i], 1/16, align=True, fontsize=9, yoffset=0.0);    i+=1
        ax.plot(intensity, bw_performance_l1_write, label="bound based on \n L1 write bandwith", color ="navy", linestyle=':')
        lines = ax.get_lines();  labelLine(lines[i], 1/16, align=True, fontsize=9, yoffset=0.0);    i+=1
        ax.plot(intensity, bw_performance_l2_rw,    label="bound based on \n L2 read+write bandwith", color ="indigo", linestyle=(0, (5, 10)))
        lines = ax.get_lines();  labelLine(lines[i], 1/40, align=True, fontsize=9, yoffset=0.0);    i+=1
        ax.plot(intensity, bw_performance_l3,       label="bound based on \n L3 read bw, write bw", color ="purple", linestyle=(0, (5, 10)))
        lines = ax.get_lines();  labelLine(lines[i], 1/40, align=True, fontsize=9, yoffset=0.0);    i+=1

    ax.plot(intensity, bw_performance_mem,      label="bound based on \n memory read+write bandwith", color ="darkviolet", linestyle='--')
    lines = ax.get_lines();  labelLine(lines[i], 1/16, align=True, fontsize=9, yoffset=0.0);    i+=1



    ax.text(1/46, 8/32,  s="x = functiondata fits into L1", fontsize="large", color="navy")
    ax.text(1/46, 5/32,  s="x = functiondata fits into L2", fontsize="large", color="indigo")
    ax.text(1/46, 3/32,  s="x = functiondata fits into L3", fontsize="large", color="purple")
    ax.text(1/46, 2/32,  s="x = functiondata doesn't fit into cache", fontsize="large", color="darkviolet")



    
    # labelLine(lines[i], 1/40, align=True, fontsize=9, yoffset=0.0) 
    # i+=1
    # labelLine(lines[i], 1/40, align=True, fontsize=9, yoffset=0.0) 
    # i+=1
    # labelLine(lines[i], 1/16, align=True, fontsize=9, yoffset=0.0) 
    # i+=1
    # labelLine(lines[i], 1/16, align=True, fontsize=9, yoffset=0.0) 
    # i+=1
    # labelLine(lines[i], 1/16, align=True, fontsize=9, yoffset=0.0)
    # i+=1

    def CPLOT(op, bt):
        return (op in op_types) and (bt in byte_types)

    # +XxoO
    for fname in fnames:
        fOI,iOI,fOIr,iOIr,fOIw,iOIw,fP,iP, nx, colors = get_data(dType, fname, L1_s, L2_s, L3_s)

        if CPLOT("flop", "both"):
            ax.plot(fOI, fP, lw=1, color=mycols[0])
            ax.scatter(fOI, fP, marker='x', s=25, lw=2, color=colors)
            ax.text( fOI[-1], 5/8*fP[-1],  s=fname+"\n flops/c,\nload+write", fontsize="xx-small", color=mycols[0])
        if CPLOT("iop", "both"):
            ax.plot(iOI, iP, lw=1, color=mycols[1])
            ax.scatter(iOI, iP, marker='x', s=25, lw=2, color=colors)
            ax.text( iOI[-1], 5/8*iP[-1],  s=fname+"\n ilops/c,\nload+write", fontsize="xx-small", color=mycols[1])
        if CPLOT("both", "both"):
            ax.plot(iOI+fOI, iP+fP, lw=1, color=mycols[2])
            ax.scatter(iOI+fOI, iP+fP, marker='x', s=25, lw=2, color=colors)
            ax.text(  (iOI+fOI)[-1], 5/8*(iP+fP)[-1],  s=fname+"\n ops/c,\nload+write", fontsize="xx-small", color=mycols[2])
        if CPLOT("flop", "load"):
            ax.plot(fOIr, fP, lw=1, color=mycols[3])
            ax.scatter(fOIr, fP, marker='x', s=25, lw=2, color=colors)
            ax.text(fOIr[-1], 5/8*fP[-1],  s=fname+"\n flops/c,\nload", fontsize="xx-small", color=mycols[3])
        if CPLOT("iop", "load"):
            ax.plot(iOIr, iP, lw=1, color=mycols[4])
            ax.scatter(iOIr, iP, marker='x', s=25, lw=2, color=colors)
            ax.text(iOIr[-1], 5/8*iP[-1],  s=fname+"\n ilops/c,\nload", fontsize="xx-small", color=mycols[4])
        if CPLOT("both", "load"):
            ax.plot(iOIr+fOIr, iP+fP, lw=1, color=mycols[5])
            ax.scatter(iOIr+fOIr,iP+fP, marker='x', s=25, lw=2, color=colors)
            ax.text((iOIr+fOIr)[-1], 5/8*(iP+fP)[-1],  s=fname+"\n ops/c,\nload", fontsize="xx-small", color=mycols[5])
        if CPLOT("flop", "store"):
            ax.plot(fOIw, fP, lw=1, color=mycols[6])
            ax.scatter(fOIw, fP, marker='x', s=25, lw=2, color=colors)
            ax.text(fOIw[-1], 5/8*fP[-1],  s=fname+"\n flops/c,\nwrite", fontsize="xx-small", color=mycols[6])
        if CPLOT("iop", "store"):
            ax.plot(iOIw, iP, lw=1, color=mycols[7])
            ax.scatter(iOIw, iP, marker='x', s=25, lw=2, color=colors)
            ax.text(iOIw[-1], 5/8*iP[-1],  s=fname+"\n ilops/c,\nwrite", fontsize="xx-small", color=mycols[7])
        if CPLOT("both", "store"):
            ax.plot(iOIw+fOIw, iP+fP, lw=1, color=mycols[8])
            ax.scatter(iOIw+fOIw, iP+fP, marker='x', s=25, lw=2, color=colors)
            ax.text((iOIw+iOIw)[-1], 5/8*(iP+fP)[-1],  s=fname+"\n ops/c,\nwrite", fontsize="xx-small", color=mycols[8])




    if(savefile):
        fig.savefig(savefile, dpi=100)
    plt.show()






def plot_kaby(dType,fnames, op_types, byte_types, mycols, sf):
        roofline(   
                    peak_performance_flops=2,
                    peak_performance_flops_fma= 4,
                    peak_performance_iops=4,
                    L1_bw=(64,32),
                    L2_bw=64,
                    L3_bw=(32,32),
                    mem_bw=16, 
                    fq=2.7,
                    L1_s=2**15, 
                    L2_s=2**18,
                    L3_s=2**24,
                    dType=dType,
                    fnames=fnames,
                    op_types=op_types,
                    byte_types=byte_types,
                    savefile = sf,
                    mycols=mycols
                )




def plot_coffee(dType,fnames, op_types, byte_types, mycols, sf):

    roofline(   
                    peak_performance_flops=2,
                    peak_performance_flops_fma= 4,
                    peak_performance_iops=4,
                    L1_bw=(64,32),
                    L2_bw=64,
                    L3_bw=(32,32),
                    mem_bw=16, 
                    fq=2.3,
                    L1_s=2**15, 
                    L2_s=2**18,
                    L3_s=2**24,
                    dType=dType,
                    fnames=fnames,
                    op_types=op_types,
                    byte_types=byte_types,
                    savefile = sf,
                    mycols=mycols
                )
                



def plot_zen3(dType,fnames, op_types, byte_types, mycols, sf):
    # 51.8gb/s -> 70?? DDR LRDDR??
    # memoryBW1 = 51.8*2**30 #B/s
    # freq = 2.7*10**9 # cycles/s
    # mm_bw = memoryBW1/freq # B/cycle = (B/s)/(cycle/s)
    # 20.6

    roofline(   
                    peak_performance_flops=4,
                    peak_performance_flops_fma= 6,
                    peak_performance_iops=4,
                    L1_bw=(32,32),
                    L2_bw=64,#(32,32),
                    L3_bw=(32,32),
                    mem_bw=20.6, 
                    fq=2.5,
                    L1_s=2**19, 
                    L2_s=2*22,
                    L3_s=2*24,
                    dType=dType,
                    fnames=fnames,
                    op_types=op_types,
                    byte_types=byte_types,
                    savefile = sf,
                    mycols=mycols,
                    zen = True
                )