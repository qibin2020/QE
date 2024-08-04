// reco analysis
Charged_Lepton get_lepton(Vec_f el_pt, Vec_f el_eta, Vec_f el_phi, Vec_i el_ch, Vec_f mu_pt, Vec_f mu_eta, Vec_f mu_phi, Vec_i mu_ch) {
    Four_Mom lep;
    int charge=0;
    // find leading electron firstly
    for (Int_t i = 0; i < el_pt.size(); i++) {
        if(el_pt.at(i)<25 || abs(el_eta.at(i))>2.5) continue;
        if(lep.Pt()<el_pt.at(i)){
            lep.SetPtEtaPhiM(el_pt.at(i),el_eta.at(i),el_phi.at(i),0.511E-3);
            charge=el_ch.at(i);
        }
    }
    // then the muon
    for (Int_t i = 0; i < mu_pt.size(); i++) {
        if(mu_pt.at(i)<25 || abs(mu_eta.at(i))>2.5) continue;
        if(lep.Pt()<mu_pt.at(i)){
            lep.SetPtEtaPhiM(mu_pt.at(i),mu_eta.at(i),mu_phi.at(i),105.65E-3);
            charge=mu_ch.at(i);
        }
    }
    return {charge,lep};
}

Four_Mom solve_nu(Four_Mom_ref lep, Vec_f nu_pt, Vec_f nu_phi) {
    Four_Mom nu;
    if(lep.Pt()<=0) return nu;
    if(nu_pt.size()==0) return nu;
    double px1,py1,pz1,px2,py2;
    px1=lep.Px();
    py1=lep.Py();
    pz1=lep.Pz();
    px2=abs(nu_pt.at(0))*cos(nu_phi.at(0));
    py2=abs(nu_pt.at(0))*sin(nu_phi.at(0));
    const double W=80.377;
    const double m=lep.M();
    double k=-((m*m) + (px1*px1) + (py1*py1) + (pz1*pz1))*(-(W*W*W*W) + 2*(W*W)*(m*m) - 4*(W*W)*px1*px2 - 4*(W*W)*py1*py2 - (m*m*m*m) + 4*(m*m)*px1*px2 + 4*(m*m)*(px2*px2) + 4*(m*m)*py1*py2 + 4*(m*m)*(py2*py2) + 4*(px1*px1)*(py2*py2) - 8*px1*px2*py1*py2 + 4*(px2*px2)*(py1*py1));
    double pz2;
    if(k>0){ // two real root
        double s1=(pz1*((W*W) - (m*m) + 2*px1*px2 + 2*py1*py2)/2 + sqrt(k)/2)/((m*m) + (px1*px1) + (py1*py1));
        double s2=(pz1*((W*W) - (m*m) + 2*px1*px2 + 2*py1*py2)/2 - sqrt(k)/2)/((m*m) + (px1*px1) + (py1*py1));
        pz2 = abs(s1) < abs(s2) ? s1 : s2;
    }else{
        pz2 = (pz1*((W*W) - (m*m) + 2*px1*px2 + 2*py1*py2)/2)/((m*m) + (px1*px1) + (py1*py1));
    }
    // build nu now
    nu.SetXYZM(px2,py2,pz2,0);
    
    return nu;
}

double get_cost_tt(Four_Mom_ref lep, Four_Mom_ref nu, Four_Mom_ref b_lep, Four_Mom_ref b_had, Four_Mom_ref j1, Four_Mom_ref j2){
    const double mt=172.76; // GeV
    const double mw=80.377;

    double c1=(lep + nu).M() - mw;
    double c2=(b_lep + lep + nu).M() - mt;
    double c3=(j1 + j2).M() - mw;
    double c4=(b_had + j1 + j2).M() - mt;

    return c1*c1 + c2*c2 +c3*c3 + c4*c4;
}

Jet_indices get_blep_bhad_j1_j2(Four_Mom_ref lep, Four_Mom_ref nu, Vec_f j_pt, Vec_f j_eta, Vec_f j_phi, Vec_f j_m, Vec_i j_btag) {
    // we find b_lep,b_had,j1,j2; follow this order
    Jet_indices ret{-1,-1,-1,-1};
    if(lep.Pt()<=0 || nu.Pt()<=0) return ret;

    // first do jet selection
    Jet_indices sel_j;
    for (Int_t i = 0; i < j_pt.size(); i++) {
        if(j_pt.at(i)<25 || abs(j_eta.at(i))>2.5) continue;
        sel_j.push_back(i);
    }
    if (sel_j.size()<4) return ret;
    // sort j by pT, note largest at the end! (to use pop)
    std::sort(sel_j.begin(), sel_j.end(), [&j_pt](const int& lhs, const int& rhs){return j_pt.at(lhs)<j_pt.at(rhs);});
    
    // 2-b jets
    Jet_indices b_j; 
    Jet_indices nonb_j; 
    for (Int_t j = 0; j < sel_j.size(); j++) { // loop from soft and hard at the end
        Int_t i=sel_j.at(j);
        if(j_btag.at(i)>0) b_j.push_back(i);
        else nonb_j.push_back(i);
    }
    // deal with differnt b-case
    // need find at least 2b, then run the matching algo
    if(b_j.size()==0){ // pop the largest one, or directly fail ?! -- ask this
        b_j.push_back(nonb_j.back());
        nonb_j.pop_back();
        b_j.push_back(nonb_j.back());
        nonb_j.pop_back();
    }else if (b_j.size()==1){ // pop the largest one
        b_j.push_back(nonb_j.back());
        nonb_j.pop_back();
    }else if (b_j.size()>2){ // keep all
    }

    // now run the matching algo.
    // do each P2, run the pseudo algo, get the mW mt mWhad mthad, minimize them
    // build permu2
    double min_cost=-1;
    Int_t ib_had_best=-1;
    Int_t ib_lep_best=-1;
    Int_t ij1_best=-1;
    Int_t ij2_best=-1;
    for(Int_t k=0;k<b_j.size();k++){
        Int_t ib_lep=b_j.at(k);
        Four_Mom b_lep;
        b_lep.SetPtEtaPhiM(j_pt.at(ib_lep),j_eta.at(ib_lep),j_phi.at(ib_lep),j_m.at(ib_lep));
        for(Int_t l=0;l<b_j.size();l++){
            Int_t ib_had=b_j.at(l);
            if(ib_had==ib_lep) continue;
            Four_Mom b_had;
            b_had.SetPtEtaPhiM(j_pt.at(ib_had),j_eta.at(ib_had),j_phi.at(ib_had),j_m.at(ib_had));
            for(Int_t m=0;m<sel_j.size()-1;m++){
                Int_t ij1=sel_j.at(m);
                if(ij1==ib_had || ij1==ib_lep) continue;
                Four_Mom j1;
                j1.SetPtEtaPhiM(j_pt.at(ij1),j_eta.at(ij1),j_phi.at(ij1),j_m.at(ij1));
                for(Int_t n=m+1;n<sel_j.size();n++){
                    Int_t ij2=sel_j.at(n);
                    if(ij2==ib_had || ij2==ib_lep) continue;
                    Four_Mom j2;
                    j2.SetPtEtaPhiM(j_pt.at(ij2),j_eta.at(ij2),j_phi.at(ij2),j_m.at(ij2));
                    double cost=get_cost_tt(lep, nu, b_lep, b_had, j1, j2);
                    if(min_cost<0 || min_cost>cost){
                        ib_had_best=ib_had;
                        ib_lep_best=ib_lep;
                        ij1_best=ij1;
                        ij2_best=ij2;
                        min_cost=cost;
                    }
                }
            }
        }
    }

    if(ib_had_best<0 || ib_lep_best<0 || ij1_best<0 || ij2_best<0) return ret;
    ret.clear();
    ret.push_back(ib_had_best);
    ret.push_back(ib_lep_best);
    ret.push_back(ij1_best);
    ret.push_back(ij2_best);
    return ret;
}

Four_Mom get_jet(int i, Vec_f j_pt, Vec_f j_eta, Vec_f j_phi, Vec_f j_m) {
    Four_Mom jet;
    if(i<0) return jet;
    jet.SetPtEtaPhiM(j_pt.at(i),j_eta.at(i),j_phi.at(i),j_m.at(i));
    return jet;
}

// now we have everything in resolved. fine.
// boosted case
Four_Mom get_Fjet(Vec_f j_pt, Vec_f j_eta, Vec_f j_phi, Vec_f j_m, Vec_f sj_pt, Vec_f sj_eta, Vec_f sj_phi, Vec_f sj_m, Vec_i sj_btag) {
    Four_Mom fjet;
    // find leading fjet
    for (Int_t i = 0; i < j_pt.size(); i++) {
        if(j_pt.at(i)<200 || abs(j_eta.at(i))>2.5) continue;
        if(fjet.Pt()<j_pt.at(i))
            fjet.SetPtEtaPhiM(j_pt.at(i),j_eta.at(i),j_phi.at(i),j_m.at(i));
    }
    return fjet;
}

bool ujet_cut(Vec_f sj_pt, Vec_f sj_eta, Vec_f sj_phi, Vec_f sj_m, Vec_i sj_btag) {
    Jet_indices sel_sj;
    for (Int_t i = 0; i < sj_pt.size(); i++) {
        if(sj_pt.at(i)<25 || abs(sj_eta.at(i))>2.5 || sj_btag.at(i)) continue;
        sel_sj.push_back(i);
    }
    if(sel_sj.size()<2) return false;
    return true;
}

double get_cost_whad(Int_t ij1, Int_t ij2, Vec_f j_pt, Vec_f j_eta, Vec_f j_phi, Vec_f j_m){
    const double mw=80.377;
    Four_Mom j1,j2;
    j1.SetPtEtaPhiM(j_pt.at(ij1),j_eta.at(ij1),j_phi.at(ij1),j_m.at(ij1));
    j2.SetPtEtaPhiM(j_pt.at(ij2),j_eta.at(ij2),j_phi.at(ij2),j_m.at(ij2));
    
    double c=(j1 + j2).M() - mw;

    return c*c;
}

double get_cost_tlep(Four_Mom_ref b_lep, Four_Mom_ref lep, Four_Mom_ref nu){
    const double mt=172.76; // GeV
    const double mw=80.377;

    double c1=(lep + nu).M() - mw; // should be fine since nu is solved
    double c2=(b_lep + lep + nu).M() - mt;

    return c1*c1 + c2*c2;
}

Jet_indices get_blep_bhad_j1_j2_boosted(Four_Mom_ref lep, Four_Mom_ref nu, Four_Mom_ref fjet, Vec_f j_pt, Vec_f j_eta, Vec_f j_phi, Vec_f j_m, Vec_i j_btag) {
    // we find b_lep,b_had,j1,j2; follow this order
    Jet_indices ret{-1,-1,-1,-1};
    if(lep.Pt()<=0 || fjet.Pt()<=0 || nu.Pt()<=0) return ret;

    // first do jet selection
    Jet_indices sel_j;
    for (Int_t i = 0; i < j_pt.size(); i++) {
        if(j_pt.at(i)<25 || abs(j_eta.at(i))>2.5) continue;
        sel_j.push_back(i);
    }
    if (sel_j.size()<4) return ret;
    // sort j by pT, note largest at the end! (to use pop)
    std::sort(sel_j.begin(), sel_j.end(), [&j_pt](const int& lhs, const int& rhs){return j_pt.at(lhs)<j_pt.at(rhs);});
    
    // had side
    // select closet 3 
    Jet_indices had_j;
    Jet_indices nonhad_j;
    for (Int_t j = 0 ; j < sel_j.size() ; j++) { // loop from soft, hardest at the last
        Int_t i=sel_j.at(j);
        Four_Mom jet;
        jet.SetPtEtaPhiM(j_pt.at(i),j_eta.at(i),j_phi.at(i),j_m.at(i));
        if(had_j.size()<3 && fjet.DeltaR(jet)<1.5)
            had_j.push_back(i); 
        else
            nonhad_j.push_back(i); 
    }
    if(had_j.size()!=3) return ret;
    // deal with b_had and j1 j2 match
    // first assum 0 1 2 and test top_had
    Int_t ij0=had_j.at(0);
    Int_t ij1=had_j.at(1);
    Int_t ij2=had_j.at(2);

    Four_Mom j0,j1,j2;
    j0.SetPtEtaPhiM(j_pt.at(ij0),j_eta.at(ij0),j_phi.at(ij0),j_m.at(ij0));
    j1.SetPtEtaPhiM(j_pt.at(ij1),j_eta.at(ij1),j_phi.at(ij1),j_m.at(ij1));
    j2.SetPtEtaPhiM(j_pt.at(ij2),j_eta.at(ij2),j_phi.at(ij2),j_m.at(ij2));
    Four_Mom top_had=j0+j1+j2;
    if(top_had.Pt()<0.9*fjet.Pt()) return ret;
    if(top_had.M()<150 || top_had.M()>225) return ret;
    Int_t nb_had=j_btag.at(ij0)+j_btag.at(ij1)+j_btag.at(ij2);
    Int_t ib_had_best=-1,ij1_best=-1,ij2_best=-1;
    if(nb_had==0){
        // select the largest, no other choice
        ib_had_best=had_j.back();
        had_j.pop_back();
        ij1_best=had_j.at(0);
        ij2_best=had_j.at(1);
    }else if(nb_had==1){
        if(j_btag.at(had_j.at(0))>0){
            ib_had_best=had_j.at(0);
            ij1_best=had_j.at(1);
            ij2_best=had_j.at(2);
        }else if (j_btag.at(had_j.at(1))>0){
            ib_had_best=had_j.at(1);
            ij1_best=had_j.at(0);
            ij2_best=had_j.at(2);
        }else{
            ib_had_best=had_j.at(2);
            ij1_best=had_j.at(0);
            ij2_best=had_j.at(1);
        }
    }else if(nb_had>=2){
        // loop b until mW_had is the best
        double c0=DBL_MAX,c1=DBL_MAX,c2=DBL_MAX; // to prevent complex logic
        if(j_btag.at(had_j.at(0))>0)
            c0=get_cost_whad(had_j.at(1),had_j.at(2),j_pt,j_eta,j_phi,j_m);
        
        if(j_btag.at(had_j.at(1))>0)
            c1=get_cost_whad(had_j.at(0),had_j.at(2),j_pt,j_eta,j_phi,j_m);

        if(j_btag.at(had_j.at(2))>0)
            c2=get_cost_whad(had_j.at(0),had_j.at(1),j_pt,j_eta,j_phi,j_m);
        // 
        if(c0<min(c1,c2)){ 
            ib_had_best=had_j.at(0); //c0 win
            ij1_best=had_j.at(1);
            ij2_best=had_j.at(2);
        }else if (c1<c2){
            ib_had_best=had_j.at(1); // c1 win
            ij1_best=had_j.at(0);
            ij2_best=had_j.at(2);
        }else{
            ib_had_best=had_j.at(2); //c2 win
            ij1_best=had_j.at(0);
            ij2_best=had_j.at(1);
        }
    }

    // deal with leptonic side
    Int_t nb_lep=0;
    Int_t ib_lep_best=-1;
    for(Int_t j=0;j<nonhad_j.size();j++){
        Int_t ib_lep=nonhad_j.at(j);
        if(j_btag.at(ib_lep)>0){
            nb_lep++;
            ib_lep_best=ib_lep; // assume the current if best
        }
    }
        
    // deal with non only 1 match
    if(nb_lep==0){
        // make leading one the b
        ib_lep_best=nonhad_j.back();
        nonhad_j.pop_back();
    }else if(nb_lep>=2){
        // ok loop nb until best match
        double min_cost=-1;
        for(Int_t j=0;j<nonhad_j.size();j++){
            Int_t ib_lep=nonhad_j.at(j);
            if(j_btag.at(ib_lep)!=1) continue;
            Four_Mom b_lep;
            b_lep.SetPtEtaPhiM(j_pt.at(ib_lep),j_eta.at(ib_lep),j_phi.at(ib_lep),j_m.at(ib_lep));
            double cost=get_cost_tlep(b_lep,lep,nu);
            if(min_cost <0 || min_cost>cost){
                ib_lep_best=ib_lep;
                min_cost=cost;
            }
        }
    }

    if(ib_lep_best <0 || ib_had_best<0 || ij1_best<0 || ij2_best<0) return ret;
    ret.clear();
    ret.push_back(ib_lep_best);
    ret.push_back(ib_had_best);
    ret.push_back(ij1_best);
    ret.push_back(ij2_best);
    
    return ret;
}      

