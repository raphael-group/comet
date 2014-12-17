
import math
import numpy as np
import scipy as sp

def tv(p, q):
    """ Total variance distance """
    return max([abs(p[i] - q[i]) for i in range(len(p))])

def discrete_convergence_eqb_plot(filelist, num_genes, ks_set, outprefix):
    
    klets_seq, gset2weight = dict(), dict()
    for f in filelist:        
        klets_seq[f] = list()
        for l in open(f+'.key'):
            v = l.rstrip().split("\t")
            for i in range(len(v)):                
                gset, w = v[i].split(":") 
                allg = set()
                for gg in gset.split(" "):
                    for g in gg.split(","):
                        allg.add(g)
                klets_seq[f].append(dict(G=",".join(sorted(allg)), W=float(w)))
                if ",".join(sorted(allg)) not in gset2weight:
                    gset2weight[",".join(sorted(allg))] = float(w)

        s_length = len(klets_seq[f])
               
    pq_list = dict()
    tv_list = dict()
    tv_w_list = dict()
    last_eqb = int(math.ceil(s_length * 0.005))
    interval_point = int(last_eqb*2) if int(last_eqb*2) > 1 else 1
    
    plot_x, plot_y, plot_yerr_low, plot_yerr_high = list(), list(), list(), list()
    y2, y2_err_h, y2_err_l = list(), list(), list()

    start_plot_index = int(s_length * 0.0) + last_eqb
    start_interval = start_plot_index % interval_point
    start_plot_index += start_interval

    print "Interval point:", interval_point
    print "Last n:", last_eqb
    print "Start plot index:", start_plot_index
        
    for i in range(start_plot_index, s_length, interval_point):
    #for i in range(last_eqb, s_length):
        
        tv_list[i] = list()
        tv_w_list[i] = list()
        #range_start = i - last_eqb + 1
        range_start = 1
        range_end = i 
        dict_f = dict()
        sum_f = dict()        
        union_f = dict()        
        union_p = list()
        union_p_w = list()

        for f in klets_seq.keys():
            #sum_f[f] -= klets_seq[f][i - last_eqb]['W']
            #sum_f[f] += klets_seq[f][i]['W']
            dict_f[f] = dict()
            sum_f[f] = dict()
            for j in range(range_start, range_end+1):               
                dict_f[f][klets_seq[f][j]['G']] = klets_seq[f][j]['W']
                if klets_seq[f][j]['G'] not in sum_f[f]:
                    #sum_f[f] += klets_seq[f][j]['W']
                    sum_f[f][klets_seq[f][j]['G']] = 0
                    #sum_w[f][klets_seq[f][j]['G']] = 0
                sum_f[f][klets_seq[f][j]['G']] += 1    
                #sum_w[f][klets_seq[f][j]['G']] += klets_seq[f][j]['W']

                if klets_seq[f][j]['G'] not in union_f:
                    union_f[klets_seq[f][j]['G']] = 0
                union_f[klets_seq[f][j]['G']] += 1
                #union_w[klets_seq[f][j]['G']] += klets_seq[f][j]['W']

        sum_union_w = sum([gset2weight[gset] for gset in union_f.keys()])
        for gset in sorted(union_f.keys()):
            union_p.append(union_f[gset]/(len(klets_seq.keys())*float(range_end-range_start+1)))        
            union_p_w.append(gset2weight[gset] / sum_union_w)


        for f in klets_seq.keys():
            #p1_dict, p2_dict = dict(), dict()            
            p1_dict = sum_f[f]
            p1_distrib = list()
            p1_distrib_w = list()
            sum_p1 = range_end - range_start + 1         
            
            sum_p1_w = sum([gset2weight[gset] if gset in sum_f[f] else 0 for gset in union_f.keys() ])
            for gset in sorted(union_f.keys()):
                if gset in sum_f[f]:
                    p1_distrib.append(sum_f[f][gset]/float(sum_p1))
                    p1_distrib_w.append(gset2weight[gset]/sum_p1_w)
                else:
                    p1_distrib.append(0)
                    p1_distrib_w.append(0)
                        
            tv_value = tv(p1_distrib, union_p)
            tv_value_w = tv(p1_distrib_w, union_p_w)        
            tv_list[i].append(tv_value)
            tv_w_list[i].append(tv_value_w)
        
        #a = mean_confidence_interval(pq_list[i])
        a2 = mean_confidence_interval(tv_list[i])
        #a = mean_confidence_interval(tv_w_list[i])
        #if i % interval_point == 0:
        plot_x.append(i)
        #plot_y.append(a[0])
        #plot_yerr_low.append(a[0] - a[1])
        #plot_yerr_high.append(a[2] - a[0])
        y2.append(a2[0])
        y2_err_h.append(a2[0] - a2[1])
        y2_err_l.append(a2[2] - a2[0])

    #plot_errorbar(plot_x, plot_y, y2, plot_yerr_low, y2_err_l, plot_yerr_high, y2_err_h, outprefix)
    #plot_errorbar(plot_x, y2, y2_err_l, y2_err_h, outprefix)

    return y2[-1]
    

def discrete_convergence(klets_seq, iter_num):
    #keys_order.append(dict(K=key, W=sum([set2scores[M]["W"] for M in collection])))
    
    tv_list = list()
    #last_eqb = int(math.ceil(s_length * 0.005))
    #interval_point = int(last_eqb*2) if int(last_eqb*2) > 1 else 1
        
    sum_num = iter_num
    
    sum_f = dict()        
    union_f = dict()        
    union_p = list()        

    for f in range(len(klets_seq)):            
        sum_f[f] = dict()
        for j in klets_seq[f].keys():                               
            sum_f[f][j] = klets_seq[f][j]['freq']
            
            if j not in union_f:
                union_f[j] = 0
            union_f[j] += klets_seq[f][j]['freq']
        
        
    for gset in sorted(union_f.keys()):
        union_p.append(union_f[gset]/(len(klets_seq)*float(sum_num)))                    

    for f in range(len(klets_seq)):        
        p1_dict = sum_f[f]
        p1_distrib = list()            

        for gset in sorted(union_f.keys()):
            if gset in sum_f[f]:
                p1_distrib.append(sum_f[f][gset]/float(sum_num))            
            else:
                p1_distrib.append(0)            
                    
        tv_value = tv(p1_distrib, union_p)
        
        tv_list.append(tv_value)
                
    a2 = mean_confidence_interval(tv_list)        
        
    return a2[0]

def discrete_convergence_check(klets_seq, s_length, conv_start):
    #keys_order.append(dict(K=key, W=sum([set2scores[M]["W"] for M in collection])))
    
    tv_list = list()
    #last_eqb = int(math.ceil(s_length * 0.005))
    #interval_point = int(last_eqb*2) if int(last_eqb*2) > 1 else 1
        
    sum_num = 0
    
    sum_f = dict()        
    union_f = dict()        
    union_p = list()        

    for f in klets_seq.keys():            
        sum_f[f] = dict()
        for j in range(len(klets_seq[f])):                               
            if klets_seq[f][j] not in sum_f[f]:
                sum_f[f][klets_seq[f][j]] = 0
        
            sum_f[f][klets_seq[f][j]] += 1    
            
            if klets_seq[f][j] not in union_f:
                union_f[klets_seq[f][j]] = 0
            union_f[klets_seq[f][j]] += 1
        sum_num = len(klets_seq[f])
                    
    for gset in sorted(union_f.keys()):
        union_p.append(union_f[gset]/(len(klets_seq.keys())*float(sum_num)))                    

    for f in klets_seq.keys():
        
        p1_dict = sum_f[f]
        p1_distrib = list()            
#        sum_p1 = range_end - range_start + 1                                 

        for gset in sorted(union_f.keys()):
            if gset in sum_f[f]:
                p1_distrib.append(sum_f[f][gset]/float(sum_num))            
            else:
                p1_distrib.append(0)            
                    
        tv_value = tv(p1_distrib, union_p)
        
        tv_list.append(tv_value)
                
    a2 = mean_confidence_interval(tv_list)        
        
    return a2[0]

#def plot_errorbar(x, y, y2, yerr_l, y2_err_l, yerr_h, y2_err_h, outprefix):
def plot_errorbar(x, y2, y2_err_l, y2_err_h, outprefix):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.figure()
    #plt.errorbar(x, y , [yerr_l, yerr_h], marker='o')
    plt.errorbar(x, y2 , [y2_err_l, y2_err_h], marker='x')
    plt.savefig(outprefix + '.freq.run.png')


def mean_confidence_interval(data, confidence=0.75):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), sp.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h
