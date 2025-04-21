import ROOT
import numpy as np
#import matplotlib.pyplot as plt
#from IPython.display import Image
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(1)
import math
import uproot
import copy
#from IPython.display import display, HTML
import time
import os

def create_legend(x1=0.15, y1=0.6, x2=0.45, y2=0.8):
    legend = ROOT.TLegend(x1, y1, x2, y2)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineWidth(0)
    legend.SetTextSize(0.05)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    return legend

def plot_histograms_with_options(histos, names, originalname, colors, options, iseff, output_name="plot.png"):
    """
    histos: dict of histograms {name: TH1}
    names: list of string names corresponding to histos
    colors: list of ROOT color codes
    options: dict of plot options (x_range, y_range, logy, ratio_plot)
    """
    #print("really plot ", histos, names, options)
    # Extract plot settings
    logy = options["logy"]
    logx = False
    if "logx" in options.keys():
        logx = options["logx"]
    ratio_plot = options["ratio_plot"]

    # Create canvas and pads if ratio plot is requested
    c = ROOT.TCanvas("c", "Canvas", 800, 800)
    if ratio_plot and len(names)>1:
        pad1 = ROOT.TPad("pad1", "Main Plot", 0, 0.3, 1, 1)
        pad1.SetBottomMargin(0.02)
        if logy:
            pad1.SetLogy()
        if logx:
            pad1.SetLogx()
        pad1.Draw()

        pad2 = ROOT.TPad("pad2", "Ratio Plot", 0, 0.0, 1, 0.31)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.3)
        pad2.Draw()
        if logx:
            pad2.SetLogx()

        pad1.cd()
    else:
        if logy:
            c.SetLogy()
        if logx:
            c.SetLogx()

    # Draw histograms
    legend = create_legend()

    maxy = -1
    miny = 999
    for i, name in enumerate(names):
        hist = histos[name]
        for binx in range(hist.GetNbinsX()):
            if hist.GetBinContent(binx+1) > maxy:
                maxy = hist.GetBinContent(binx+1)
            if hist.GetBinContent(binx+1) < miny and hist.GetBinContent(binx+1)>0 :
                miny = hist.GetBinContent(binx+1)

    for i, name in enumerate(names):
        hist = histos[name]
        print(hist.Integral())
        hist.SetLineColor(colors[i])
        hist.SetLineWidth(2)

        hist.GetYaxis().SetNdivisions(505)
        hist.GetYaxis().SetTitleSize(30)
        hist.GetYaxis().SetTitleFont(43)
        hist.GetYaxis().SetTitleOffset(1.4)
        hist.GetYaxis().SetLabelFont(43)
        hist.GetYaxis().SetLabelSize(25)

        if not(ratio_plot and len(names)>1):
            hist.GetXaxis().SetTitleSize(30)
            hist.GetXaxis().SetTitleFont(43)
            hist.GetXaxis().SetTitleOffset(1.0)
            hist.GetXaxis().SetLabelFont(43)
            hist.GetXaxis().SetLabelSize(25)

        # Axis ranges
        hist.GetXaxis().SetRangeUser(options['x_range'][0],options['x_range'][1])
        hist.GetYaxis().SetRangeUser(options['y_range'][0],options['y_range'][1])
        if logy:
            hist.GetYaxis().SetRangeUser(miny*0.1, maxy*10)
        else:
            hist.GetYaxis().SetRangeUser(miny*0.8, maxy*1.2)
        hist.GetXaxis().SetTitle(options["xlabel"])
        hist.GetYaxis().SetTitle(options["ylabel"])
        #print(options['x_range'], options['y_range'])
        #print(options['xlabel'], options['ylabel'])
        if "DIST2D" in hist.GetName():
            draw_opt = "colz" if i == 0 else "colz SAME"
        else:
            draw_opt = "HIST" if i == 0 else "HIST SAME"
        hist.Draw(draw_opt)
        legend.AddEntry(hist, originalname[i], "L")

    legend.Draw()

    # Ratio plot
    #print(names)
    if ratio_plot and len(names)>1:
        pad2.cd()

        # Take first histogram as reference
        h_ref = histos[names[0]].Clone("ref")
        h_ref.SetLineColor(ROOT.kBlack)
        h_ref.SetLineWidth(2)

        # Ratio histograms
        for i, name in enumerate(names):
            if i==0: continue
            if iseff:
                ratio = histos[name].Clone("ratio_"+name)
                ratio.Divide(h_ref)
            else:
                ratio = h_ref.Clone("ratio_"+name)
                ratio2 = histos[name].Clone("ratio_"+name)
                ratio2.Add(h_ref)
                ratio.Divide(ratio2)
            ratio.SetLineColor(colors[i])
            ratio.SetMarkerStyle(22)
            ratio.SetMarkerColor(colors[i])
            ratio.SetLineWidth(2)
            if i == 1:
                ratio.GetYaxis().SetNdivisions(505)
                ratio.GetYaxis().SetTitleSize(30)
                ratio.GetYaxis().SetTitleFont(43)
                ratio.GetYaxis().SetTitleOffset(1.4)
                ratio.GetYaxis().SetLabelFont(43)
                ratio.GetYaxis().SetLabelSize(25)

                ratio.GetXaxis().SetTitleSize(30)
                ratio.GetXaxis().SetTitleFont(43)
                ratio.GetXaxis().SetTitleOffset(1.0)
                ratio.GetXaxis().SetLabelFont(43)
                ratio.GetXaxis().SetLabelSize(25)

                ratio.GetXaxis().SetRangeUser(options['x_range'][0],options['x_range'][1])
                ratio.GetYaxis().SetRangeUser(options['ratio_ylim'][0],options['ratio_ylim'][1])
                ratio.GetXaxis().SetTitle(options["xlabel"])
                ratio.GetYaxis().SetTitle("Eff.")
                ratio.Draw("P")
            else:
                ratio.GetYaxis().SetNdivisions(505)
                ratio.GetYaxis().SetTitleSize(30)
                ratio.GetYaxis().SetTitleFont(43)
                ratio.GetYaxis().SetTitleOffset(1.4)
                ratio.GetYaxis().SetLabelFont(43)
                ratio.GetYaxis().SetLabelSize(25)

                ratio.GetXaxis().SetTitleSize(30)
                ratio.GetXaxis().SetTitleFont(43)
                ratio.GetXaxis().SetTitleOffset(1.0)
                ratio.GetXaxis().SetLabelFont(43)
                ratio.GetXaxis().SetLabelSize(25)

                ratio.GetXaxis().SetRangeUser(options['x_range'][0],options['x_range'][1])
                ratio.GetYaxis().SetRangeUser(options['ratio_ylim'][0],options['ratio_ylim'][1])
                ratio.GetXaxis().SetTitle(options["xlabel"])
                ratio.GetYaxis().SetTitle("Eff.")
                ratio.Draw("P SAME")
    c.SaveAs(output_name+'.png')
    c.SaveAs(output_name+'.C')

def create_histograms(plot_options, name_suffix):
    """Create and return a dictionary of histograms with the given suffix."""
    rename = name_suffix.replace(" ", "").replace("(", "").replace(")", "").replace(",", "")
    histos = {}
    for key in plot_options:
        if "DIST2D" in key:
            histos[key] = ROOT.TH2F(key+"_"+rename, ';'+plot_options[key]["xlabel"]+';'+plot_options[key]["ylabel"], len(plot_options[key]["bin_edges"]) - 1, plot_options[key]["bin_edges"], len(plot_options[key]["ybin_edges"]) - 1, plot_options[key]["ybin_edges"])
        else:
            histos[key] = ROOT.TH1F(key+"_"+rename, ';'+plot_options[key]["xlabel"]+';'+plot_options[key]["ylabel"], len(plot_options[key]["bin_edges"]) - 1, plot_options[key]["bin_edges"])
    return histos

def normalize_histogram(hist, min_y=0.01):
    """Scale a histogram to unit area and set drawing properties."""
    if hist.Integral() != 0:
        hist.Scale(1.0 / hist.Integral())
    hist.GetYaxis().SetRangeUser(min_y, 1)


def delta_phi(phi1, phi2):
    """
    Compute the difference in phi between two angles,
    wrapping it into the interval [-pi, pi].
    """
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi <= -math.pi:
        dphi += 2 * math.pi
    return dphi

def delta_r(eta1, phi1, eta2, phi2):
    """
    Compute delta R between two objects with eta and phi.
    """
    d_eta = eta1 - eta2
    d_phi = delta_phi(phi1, phi2)
    return math.sqrt(d_eta**2 + d_phi**2)

def delta_r_numpy(eta1, phi1, eta2_array, phi2_array):
    """Compute the delta R between a single (eta1, phi1) and multiple (eta2_array, phi2_array)."""
    delta_eta = eta1 - eta2_array
    delta_phi = np.arctan2(np.sin(phi1 - phi2_array), np.cos(phi1 - phi2_array))  # 
    return np.sqrt(delta_eta**2 + delta_phi**2)

# B
         
def fill_histograms(root_file, branchtoread, TOBEFILL, histos, plot_options, rename, cuts, isgen):
    with uproot.open(root_file) as file:
        tree = file["Events"]  # Replace with the actual tree name
        # Load all relevant branches at once
        data = tree.arrays(branchtoread, library="np")

    # Filter events with nKshort > 0
    valid_events = len(data["nB"])
    for ievt in range(valid_events):
        genmatch_bpion_idx = data["Genmatch_idx_bpion"][ievt]
        genmatch_b_idx = data["Genmatch_idx_b"][ievt]
        genmatch_D0_idx = data["Genmatch_idx_d0"][ievt]
        genmatch_ks0_idx = data["Genmatch_idx_ks"][ievt]
        gen_eta = data["GenPart_eta"][ievt]
        gen_phi = data["GenPart_phi"][ievt]
        gen_pt = data["GenPart_pt"][ievt]
        B_Kin_pt       =    data["B_B_Kin_pt"][ievt]
        B_Kin_eta      =    data["B_B_Kin_eta"][ievt]
        B_Kin_phi      =    data["B_B_Kin_phi"][ievt]
        D_Kin_pt       =    data["B_D0_Kin_pt"][ievt]
        D_Kin_eta      =    data["B_D0_Kin_eta"][ievt]
        D_Kin_phi      =    data["B_D0_Kin_phi"][ievt]
        closest_dR = []
        closest_idx = []
        idx_d0 = genmatch_D0_idx[0]
        if idx_d0<0: continue
        idx_B = genmatch_b_idx[0]
        if idx_B<0: continue
        if len(D_Kin_pt)<=0: continue

        # match B
        gen_eta_i, gen_phi_i = gen_eta[idx_B], gen_phi[idx_B]
        dR_values = delta_r_numpy(gen_eta_i, gen_phi_i, B_Kin_eta, B_Kin_phi)
        # match D
        #gen_eta_i, gen_phi_i = gen_eta[idx_d0], gen_phi[idx_d0]
        #dR_values = delta_r_numpy(gen_eta_i, gen_phi_i, D_Kin_eta, D_Kin_phi)

        sorted_indices = np.argsort(dR_values)[:2]  # Get two smallest dR
        closest_dR.append(dR_values[sorted_indices])
        closest_idx.append(sorted_indices)
        idx_reco = closest_idx[0][0]
        sorted_indices_others = np.argsort(dR_values)[1:]

        delta_pt = np.abs(B_Kin_pt[idx_reco] - gen_pt[idx_B]) / gen_pt[idx_B]
        #delta_pt = np.abs(D_Kin_pt[idx_reco] - gen_pt[idx_d0]) / gen_pt[idx_d0]
        dR_cut = 0.1
        dpt_cut = 0.2    
        histos["DIST2D_dR_dpt"].Fill(abs(delta_pt), closest_dR[0][0])
        #print(abs(delta_pt), closest_dR[0][0])
        if (( cuts == "matched" and closest_dR[0][0]<dR_cut and abs(delta_pt)<dpt_cut )
            or (cuts == "unmatched" and not(abs(delta_pt)<dpt_cut and closest_dR[0][0]<dR_cut))
            or ("comb.bkg" in cuts and closest_dR[0][0]<dR_cut and abs(delta_pt)<dpt_cut)
           ):            
            if "comb.bkg" in cuts:
                for it in TOBEFILL:
                    if it=="DIST2D_dR_dpt":
                        continue
                    for j in range(len(sorted_indices_others)):
                        idx_reco_other = sorted_indices_others[j]
                        if "_DD_" in it:
                            name1=it.split("_DD_")[0]
                            name2=it.split("_DD_")[1]
                            if data[name2][ievt][idx_reco_other]!=0:
                                histos[it].Fill( data[name1][ievt][idx_reco_other]/data[name2][ievt][idx_reco_other] )
                        elif "_QQ_" in it:
                            name1=it.split("_QQ_")[0]
                            name2=it.split("_QQ_")[1]
                            histos[it].Fill( np.sqrt(data[name1][ievt][idx_reco_other]**2+data[name2][ievt][idx_reco_other]**2) )
                        elif "Kshort_sigmaDistMagXY" in it:
                            histos[it].Fill(1.0/data[it][ievt][idx_reco_other])
                        else:
                            histos[it].Fill(data[it][ievt][idx_reco_other])
                
            else:
                for it in TOBEFILL:
                    if it=="DIST2D_dR_dpt":
                        continue
                    if "_DD_" in it:
                        name1=it.split("_DD_")[0]
                        name2=it.split("_DD_")[1]
                        if data[name2][ievt][idx_reco]!=0:
                            histos[it].Fill( data[name1][ievt][idx_reco]/data[name2][ievt][idx_reco] )
                    elif "_QQ_" in it:
                        name1=it.split("_QQ_")[0]
                        name2=it.split("_QQ_")[1]
                        histos[it].Fill( np.sqrt(data[name1][ievt][idx_reco]**2+data[name2][ievt][idx_reco]**2) )
                    elif "Kshort_sigmaDistMagXY" in it:
                        histos[it].Fill(1.0/data[it][ievt][idx_reco])
                    else:
                        histos[it].Fill(data[it][ievt][idx_reco])

def define_hist_setting(TOBEFILL):
    plot_option_gen = {}                            
    binlist={}
    rangelist={}
    ifratio_plot=True
    for it in TOBEFILL:
        default_logx=False
        default_logy=True
        default_xlabel=it
        if "_pt" in it:
            binlist[it]=(0.01, 50)
            #rangelist[it]=np.linspace(0, 15, 100)
            rangelist[it] = np.logspace(np.log10(0.01), np.log10(50), 100)
            default_logx=True
            default_logy=True
            default_xlabel=it+" [GeV]"
        elif "_eta" in it:
            binlist[it]=(-4, 4)
            rangelist[it]=np.linspace(-4, 4, 100)
        elif "_phi" in it:
            binlist[it]=(-3.2, 3.2)
            rangelist[it]=np.linspace(-3.2, 3.2, 100)
        elif "Kshort_massSquared" == it:
            binlist[it]=(0.15, 0.36)
            rangelist[it]=np.linspace(0.15, 0.36, 100)
            default_xlabel=it+" [GeV]"
        elif "D0_mass" in it or "D0_Kin_mass" in it or "B_D0_Kinfitmass" in it:
            binlist[it]=(1.86484-0.3, 1.86484+0.3)
            rangelist[it]=np.linspace(1.86484-0.3, 1.86484+0.3, 100)
            default_xlabel=it+" [GeV]"
        elif "B_premass" in it or "B_Kinfitmass" in it or "B_Kin_mass" in it:
            binlist[it]=(5.27934-0.5, 5.27934+0.5)
            rangelist[it]=np.linspace(5.27934-0.5, 5.27934+0.5, 100)
            default_xlabel=it+" [GeV]"
        elif "_mass" in it:
            binlist[it]=(0.35, 0.62)
            rangelist[it]=np.linspace(0.35, 0.62, 100)
            default_xlabel=it+" [GeV]"
        elif "cxPtR2" in it:
            binlist[it]=(0.001, 10000)
            #rangelist[it] = np.linspace(0, 10000, 100)
            rangelist[it] = np.logspace(np.log10(0.001), np.log10(10000), 100)
            default_logx=True
            default_logy=True
            default_xlabel=it+" [cm]"
        elif "cxPtz" in it:
            binlist[it]=(0.01, 100)
            #rangelist[it]=np.linspace(0, 100, 100)
            rangelist[it] = np.logspace(np.log10(0.01), np.log10(100), 100)
            default_logx=True
            default_logy=True
            default_xlabel=it+" [cm]"
        elif "_dca" in it:
            binlist[it]=(0.0001, 100)
            #rangelist[it]=np.linspace(0, 1, 100)
            rangelist[it] = np.logspace(np.log10(0.0001), np.log10(100), 100)
            default_xlabel=it+" [cm]"
            default_logx=True
        elif "Kshort_sigmaDistMagXYZ" in it:
            binlist[it]=(0, 10)
            rangelist[it]=np.linspace(0, 10, 100)
            default_xlabel="Ks0 #sigma(l_{xyz})"
        elif "Kshort_sigmaDistMagXY" in it:
            binlist[it]=(0, 10)
            rangelist[it]=np.linspace(0, 10, 100)
            default_xlabel="Ks0 #sigma(l_{xy})"
        elif "D0_Kin_l_xy" in it:
            binlist[it]=(0, 10)
            rangelist[it]=np.linspace(0, 10, 100)
            default_xlabel="D0 #sigma(l_{xy})"        
        elif "Kshort_trk1_dot_trk2" in it or "D0_trk3_dot_trk4" in it or "dot" in it:
            binlist[it]=(0.01, 100)
            #rangelist[it]=np.linspace(-1, 20, 100)
            rangelist[it] = np.logspace(np.log10(0.01), np.log10(100), 100)
            default_logx=True
        elif "_chi2" in it or "Kshort_KLM_normalizedChi2" in it:
            binlist[it]=(-1, 15)
            rangelist[it]=np.linspace(-1, 15, 110)
        elif "_prob" in it:
            binlist[it]=(0, 1)
            rangelist[it]=np.linspace(0, 1, 110)
        elif "Ks0_Kin_dof" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
        elif "dof" in it:
            binlist[it]=(0, 100)
            rangelist[it]=np.linspace(0, 100, 100)
        elif "cos_theta" in it:
            #binlist[it]=(0, 1)
            #rangelist[it]=np.linspace(0, 1, 110)
            binlist[it]=(0.9, 1)
            rangelist[it] = np.logspace(np.log10(0.9), np.log10(1), 100)
            default_logx=True
        elif "_DD_" in it:
            binlist[it]=(0.1, 100)
            #rangelist[it]=np.linspace(0, 100, 100)
            rangelist[it] = np.logspace(np.log10(0.1), np.log10(100), 100)
            if "Kshort_Kin_l_xy_DD_Kshort_Kin_l_xy_unc" in it:
                default_xlabel="Ks0 KinFit #sigma(l_{xy})"
            if "Kshort_KLM_distMagXY_DD_Kshort_KLM_sigmaDistMagXY" in it:
                default_xlabel="Ks0 KLMFit #sigma(l_{xy})"
            if "Kshort_KLM_distMagXYZ_DD_Kshort_KLM_sigmaDistMagXYZ" in it:
                default_xlabel="Ks0 KLMFit #sigma(l_{xyz})"
            default_logx=True
        elif "vtx_y" in it:
            binlist[it]=(0.01, 50)
            #rangelist[it]=np.linspace(0, 50, 100)
            rangelist[it] = np.logspace(np.log10(0.01), np.log10(50), 100)
            default_logx=True
            default_xlabel=it.split("_QQ_")[0].replace("_x","_xy")+" [cm]"
        elif "vtx_z" in it:
            binlist[it]=(0.01, 50)
            #rangelist[it]=np.linspace(0, 50, 100)
            rangelist[it] = np.logspace(np.log10(0.01), np.log10(50), 100)
            default_logx=True
            default_xlabel=it+" [cm]"
        elif "flight_distance" in it:
            binlist[it]=(0.01, 100)
            rangelist[it] = np.logspace(np.log10(0.01), np.log10(100), 100)
            default_logx=True
            default_xlabel=it+" [cm]"
        elif "Trk1" in it and "ipsigXY" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel="Track1 ipsigXY"
        elif "Trk2" in it and "ipsigXY" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel="Track2 ipsigXY"
        elif "Trk3" in it and "ipsigXY" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel="Track3 ipsigXY"
        elif "Trk4" in it and "ipsigXY" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel="Track4 ipsigXY"
        elif "Trk1" in it and "ipsigZ" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel="Track1 ipsigZ"
        elif "Trk2" in it and "ipsigZ" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel="Track2 ipsigZ"
        elif "Trk3" in it and "ipsigZ" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel="Track3 ipsigZ"
        elif "Trk4" in it and "ipsigZ" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel="Track4 ipsigZ"
        elif "normalizedChi2" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel=it
        elif "TrackHits" in it:
            binlist[it]=(0, 50)
            rangelist[it]=np.linspace(0, 50, 100)
            default_xlabel=it
        elif "PixelHits" in it:
            binlist[it]=(0, 20)
            rangelist[it]=np.linspace(0, 20, 100)
            default_xlabel=it
        if "DIST2D_dR_dpt" == it: continue 
        plot_option_gen[it] =  {
            "x_range": binlist[it],
            "y_range": (1e-2, 5),
            "ratio_ylim": (0, 1.5),
            "logx": default_logx,
            "logy": default_logy,
            "ratio_plot": ifratio_plot,
            'xvar': "", 
            "xlabel": default_xlabel,
            "ylabel": "A.U.",
            "bin_edges": rangelist[it],
        }

    plot_option_gen["DIST2D_dR_dpt"] = {
            "x_range": (0, 1),
            "y_range": (0, 3.0),
            "ratio_ylim": (0, 1.5),
            "logx": False,
            "logy": False,
            "ratio_plot": False,
            'xvar': "",
            "xlabel": "delta(pT)",
            "ylabel": "delta(R)",
            "bin_edges": np.linspace(0, 1, 50),
            "ybin_edges": np.linspace(0, 3.0, 50)
    }
    return plot_option_gen

def main():

    filename_list=[
        #"/eos/uscms/store/user/yilai/postprocess_BPHNano1.root",
        #"/eos/uscms/store/user/yilai/postprocess_BPHNano1.root",
        "/eos/uscms/store/user/yilai/postprocess_BPHNano_all.root",
        "/eos/uscms/store/user/yilai/postprocess_BPHNano_all.root",
        #"test.root",
        #"test.root"
    ]
    
    names=[
        "matched",
        "matched (comb.bkg)",
    ]
    
    cuts = [
            "matched",
            "matched (comb.bkg)",
           ]

    TOBEFILL_set1=[
                   "B_B_premass", "B_B_Kinfitmass",
                   "B_B_Kin_vtx_x_QQ_B_B_Kin_vtx_y", "B_B_Kin_vtx_z", "B_B_Kin_chi2", "B_B_Kin_dof", "B_B_Kin_prob",
                   "B_B_Kin_mass",
                   "B_B_Kin_fitted_cos_theta_2D", "B_B_Kin_l_xy_DD_B_B_Kin_l_xy_unc",
                   "B_B_Kin_trk_pt", "B_B_Kin_trk_eta", "B_B_Kin_D0_pt", "B_B_Kin_D0_eta",
                   "B_D0_flight_distance", "B_D0_flight_distance_2D",
                   ]
    TOBEFILL_set2=[
                   "B_Ks0_Kin_vtx_x_QQ_B_Ks0_Kin_vtx_y", "B_Ks0_Kin_vtx_z", "B_Ks0_Kin_chi2", "B_Ks0_Kin_dof", "B_Ks0_Kin_prob",
                   "B_Ks0_Kin_pt", "B_Ks0_Kin_eta", "B_Ks0_Kin_mass", "B_Ks0_Kin_fitted_cos_theta_2D", "B_Ks0_Kin_l_xy_DD_B_Ks0_Kin_l_xy_unc",
                   "B_Ks0_Kin_trk1_pt", "B_Ks0_Kin_trk1_eta", "B_Ks0_Kin_trk2_pt", "B_Ks0_Kin_trk2_eta",
                   ] 
    TOBEFILL_set3=[
                   "B_D0_Kinfitmass", "B_D0_Kin_vtx_x_QQ_B_D0_Kin_vtx_y", "B_D0_Kin_vtx_z", "B_D0_Kin_chi2", "B_D0_Kin_dof", "B_D0_Kin_prob",
    "B_D0_Kin_mass", "B_D0_Kin_fitted_cos_theta_2D", "B_D0_Kin_l_xy_DD_B_D0_Kin_l_xy_unc",
                   "B_D0_Kin_trk3_pt", "B_D0_Kin_trk3_eta", "B_D0_Kin_trk4_pt", "B_D0_Kin_trk4_eta",
                   "B_D0_Kin_ks0_pt", "B_D0_Kin_ks0_eta", "B_ks0_flight_distance", "B_ks0_flight_distance_2D",
             ]
    TOBEFILL_set4=[
            "B_DiTrk1_cxPtR2", "B_DiTrk1_cxPtz", "B_DiTrk1_dot", "B_DiTrk1_trk1_dcaSig", "B_DiTrk1_trk2_dcaSig",
            "B_DiTrk1_dca", "B_DiTrk1_massSquared", "B_DiTrk1_KLM_vtx_x_QQ_B_DiTrk1_KLM_vtx_y", "B_DiTrk1_KLM_vtx_z",
            ]
    TOBEFILL_set5=[
            "B_DiTrk1_KLM_chi2", "B_DiTrk1_KLM_ndof", "B_DiTrk1_KLM_normalizedChi2", "B_DiTrk1_KLM_distMagXY_DD_B_DiTrk1_KLM_distMagXYErr", "B_DiTrk1_KLM_cos_theta_XY",
            ]
    TOBEFILL_set6=[
                "B_DiTrk2_cxPtR2", "B_DiTrk2_cxPtz", "B_DiTrk2_dot", "B_DiTrk2_trk2_dcaSig", "B_DiTrk2_trk2_dcaSig",
            "B_DiTrk2_dca", "B_DiTrk2_massSquared", "B_DiTrk2_KLM_vtx_x_QQ_B_DiTrk2_KLM_vtx_y", "B_DiTrk2_KLM_vtx_z"]
    TOBEFILL_set7=[
            "B_DiTrk2_KLM_chi2", "B_DiTrk2_KLM_ndof", "B_DiTrk2_KLM_normalizedChi2", "B_DiTrk2_KLM_distMagXY_DD_B_DiTrk2_KLM_distMagXYErr", "B_DiTrk2_KLM_cos_theta_XY",
            ]

    branch2read_set1 = [
                   "B_B_premass", "B_B_Kinfitmass",
                   "B_B_Kin_vtx_x", "B_B_Kin_vtx_y", "B_B_Kin_vtx_z", "B_B_Kin_chi2", "B_B_Kin_dof", "B_B_Kin_prob",
                   "B_B_Kin_mass",
                   "B_B_Kin_fitted_cos_theta_2D", "B_B_Kin_l_xy", "B_B_Kin_l_xy_unc",
                   "B_B_Kin_trk_pt", "B_B_Kin_trk_eta", "B_B_Kin_D0_pt", "B_B_Kin_D0_eta",
                   "B_D0_flight_distance", "B_D0_flight_distance_2D",
                   ]
    branch2read_set2 = [
                   "B_Ks0_Kin_vtx_x", "B_Ks0_Kin_vtx_y", "B_Ks0_Kin_vtx_z", "B_Ks0_Kin_chi2", "B_Ks0_Kin_dof", "B_Ks0_Kin_prob",
                   "B_Ks0_Kin_pt", "B_Ks0_Kin_eta", "B_Ks0_Kin_mass", "B_Ks0_Kin_fitted_cos_theta_2D", "B_Ks0_Kin_l_xy", "B_Ks0_Kin_l_xy_unc",
                   "B_Ks0_Kin_trk1_pt", "B_Ks0_Kin_trk1_eta", "B_Ks0_Kin_trk2_pt", "B_Ks0_Kin_trk2_eta",
                   ]
    branch2read_set3 = [
                   "B_D0_Kinfitmass", "B_D0_Kin_vtx_x", "B_D0_Kin_vtx_y", "B_D0_Kin_vtx_z", "B_D0_Kin_chi2", "B_D0_Kin_dof", "B_D0_Kin_prob",
                   "B_D0_Kin_mass", "B_D0_Kin_fitted_cos_theta_2D", "B_D0_Kin_l_xy", "B_D0_Kin_l_xy_unc",
                   "B_D0_Kin_trk3_pt", "B_D0_Kin_trk3_eta", "B_D0_Kin_trk4_pt", "B_D0_Kin_trk4_eta",
                   "B_D0_Kin_ks0_pt", "B_D0_Kin_ks0_eta", "B_ks0_flight_distance", "B_ks0_flight_distance_2D",
                   ]
    branch2read_set4 = [
                        "B_DiTrk1_cxPtR2", "B_DiTrk1_cxPtz", "B_DiTrk1_dot", "B_DiTrk1_trk1_dcaSig", "B_DiTrk1_trk2_dcaSig",
            "B_DiTrk1_dca", "B_DiTrk1_massSquared", "B_DiTrk1_KLM_vtx_x", "B_DiTrk1_KLM_vtx_y", "B_DiTrk1_KLM_vtx_z"
    ]
    branch2read_set5 = [
            "B_DiTrk1_KLM_chi2", "B_DiTrk1_KLM_ndof", "B_DiTrk1_KLM_normalizedChi2", "B_DiTrk1_KLM_distMagXY", "B_DiTrk1_KLM_distMagXYErr", "B_DiTrk1_KLM_cos_theta_XY",
            ]
    branch2read_set6 = [
            "B_DiTrk2_cxPtR2", "B_DiTrk2_cxPtz", "B_DiTrk2_dot", "B_DiTrk2_trk2_dcaSig", "B_DiTrk2_trk2_dcaSig",
            "B_DiTrk2_dca", "B_DiTrk2_massSquared", "B_DiTrk2_KLM_vtx_x", "B_DiTrk2_KLM_vtx_y", "B_DiTrk2_KLM_vtx_z"]
    branch2read_set7 = [
            "B_DiTrk2_KLM_chi2", "B_DiTrk2_KLM_ndof", "B_DiTrk2_KLM_normalizedChi2", "B_DiTrk2_KLM_distMagXY", "B_DiTrk2_KLM_distMagXYErr", "B_DiTrk2_KLM_cos_theta_XY",

            ]

    TOBEFILL_must_have = ["DIST2D_dR_dpt" ]
    branch2read_must_have = [
                   "GenPart_eta", "GenPart_phi", "GenPart_pt",
                   "nGenmatch", "Genmatch_idx_bpion", "Genmatch_idx_b", "Genmatch_idx_d0", "Genmatch_idx_ks",
                   "nB",
                   "B_B_Kin_pt", "B_B_Kin_eta", "B_B_Kin_phi",
                   "B_D0_Kin_pt", "B_D0_Kin_eta", "B_D0_Kin_phi"]

    TOBEFILL = TOBEFILL_must_have + TOBEFILL_set2
    branch2read = branch2read_must_have + branch2read_set2
    TOBEFILL = TOBEFILL_must_have + TOBEFILL_set8
    branch2read = branch2read_must_have + branch2read_set8


    plot_options = define_hist_setting(TOBEFILL)

    # Colors for each plot
    colors = [
        ROOT.kRed,
        ROOT.kBlue,
        ROOT.kOrange,
        ROOT.kCyan,
        ROOT.kViolet,
        ROOT.kBlack,
    ]

    # Load files and trees
    rawfile_list = [fname for fname in filename_list]
    file_list = [ROOT.TFile.Open(fname) for fname in filename_list]
    tree_list = [f.Get("Events") for f in file_list]

    all_histograms = {}

    renamelist=[]
    for index, name in enumerate(names):
        print(index, name)
        rename = name.replace(" ", "").replace("(", "").replace(")", "").replace(",", "")
        renamelist.append(rename)
        histos = create_histograms(plot_options, rename)
        fill_histograms(rawfile_list[index], branch2read, TOBEFILL, histos, plot_options, rename, cuts[index], True)
        for h in histos.keys():
            normalize_histogram(histos[h])
        all_histograms[rename] = histos

    os.system("mkdir -p pic_recoB")
    # Plot each variable with its options
    for var, opts in plot_options.items():
        histos_to_plot = {name: all_histograms[name][var] for name in renamelist}
        iseff=True
        plot_histograms_with_options(histos_to_plot, renamelist, names, colors, opts, iseff, output_name=f"pic_recoB/{var}")


main()

