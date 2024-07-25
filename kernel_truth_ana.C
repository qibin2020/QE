int find_particle(int psize, int target_pid, int target_status, Vec_i PID, Vec_i STATUS){
    for (Int_t i = 0; i < psize; i++) {
        if((PID.at(i)==target_pid) && STATUS.at(i)==target_status)
            return i;
    }
    return -1; // not possible
}

int get_top(int psize, bool anti, Vec_i PID, Vec_i STATUS){
    return find_particle(psize,anti?-6:6,22,PID,STATUS);
}

int get_W(int psize, bool neg, Vec_i PID, Vec_i STATUS){
    return find_particle(psize,neg?-24:24,22,PID,STATUS);
}

int get_b(int psize, bool neg, Vec_i PID, Vec_i STATUS){
    return find_particle(psize,neg?-5:5,23,PID,STATUS);
}

int find_lep(int psize, Vec_i PID, Vec_i STATUS, Vec_i M1, Vec_i M2){
    for (Int_t i = 0; i < psize; i++) {
        if((abs(PID.at(i))==11 || abs(PID.at(i))==13) 
                && (STATUS.at(i)==1 || STATUS.at(i)==23)
                && ((M1.at(i)>=0 && abs(PID.at(M1.at(i)))==24) || (M2.at(i)>=0 && abs(PID.at(M2.at(i)))==24)))
            return i;
    }
    return -1; // not possible
}

int find_nu(int psize, Vec_i PID, Vec_i STATUS, Vec_i M1, Vec_i M2){
    for (Int_t i = 0; i < psize; i++) {
        if((abs(PID.at(i))==12 || abs(PID.at(i))==14) 
                && (STATUS.at(i)==1 || STATUS.at(i)==23)
                && ((M1.at(i)>=0 && abs(PID.at(M1.at(i)))==24) || (M2.at(i)>=0 && abs(PID.at(M2.at(i)))==24)))
            return i;
    }
    return -1; // not possible
}

int find_u(int psize, Vec_i PID, Vec_i STATUS, Vec_i M1, Vec_i M2){
    for (Int_t i = 0; i < psize; i++) {
        if((abs(PID.at(i))==2 || abs(PID.at(i))==4) 
                && STATUS.at(i)==23
                && ((M1.at(i)>=0 && abs(PID.at(M1.at(i)))==24) || (M2.at(i)>=0 && abs(PID.at(M2.at(i)))==24)))
            return i;
    }
    return -1; // not possible
}

int find_d(int psize, Vec_i PID, Vec_i STATUS, Vec_i M1, Vec_i M2){
    for (Int_t i = 0; i < psize; i++) {
        if((abs(PID.at(i))==1 || abs(PID.at(i))==3) 
                && STATUS.at(i)==23
                && ((M1.at(i)>=0 && abs(PID.at(M1.at(i)))==24) || (M2.at(i)>=0 && abs(PID.at(M2.at(i)))==24)))
            return i;
    }
    return -1; // not possible
}


Four_Mom get_particle_idx(int i, Vec_d E, Vec_d PT, Vec_d Eta, Vec_d Phi){
    Four_Mom par;
    par.SetPtEtaPhiE(PT.at(i), Eta.at(i), Phi.at(i), E.at(i));
    return par;
}


// for debug
double get_theta_jet(std::string axis, Four_Mom_ref _had_top, Frame_ref helix, Four_Mom_ref ttbar, Four_Mom_ref par){ 
    Four_Mom had_top = _had_top;
    had_top.Boost(-ttbar.BoostVector());
    
    Four_Mom jet = par;
    jet.Boost(-ttbar.BoostVector());
    jet.Boost(-had_top.BoostVector());

    return helix.at(axis).Dot(jet.Vect().Unit());
}

