Frame get_top_frame(Four_Mom_ref _top, Four_Mom_ref ttbar) {
    //ok, first boost in ttbar CM
    Four_Mom top = _top;
    top.Boost(-ttbar.BoostVector());
    
    // Calculate helicity basis vectors
    Three_Mom k_hat = top.Vect().Unit();
    
    // Define beam direction
    Three_Mom p_hat(0, 0, 1);
    //Three_Mom p_hat=ttbar.Vect().Unit(); // beam axis should be parallel to z!
    
    double y_p = p_hat.Dot(k_hat);
    double r_p = TMath::Sqrt(1 - y_p*y_p);  
    
    // Calculate r_hat
    Three_Mom r_hat = 1 / r_p * (p_hat - y_p * k_hat);
    r_hat = r_hat.Unit();
    
    // Calculate n_hat
    Three_Mom n_hat = r_hat.Cross(k_hat);
    n_hat = n_hat.Unit();

    return {{"k", k_hat}, {"r", r_hat}, {"n", n_hat}};
}

double get_cos_theta(Four_Mom_ref _top, Four_Mom_ref ttbar) {
    //ok, first boost in ttbar CM
    Four_Mom top = _top;
    top.Boost(-ttbar.BoostVector());
    
    // Calculate helicity basis vectors
    Three_Mom k_hat = top.Vect().Unit();
    
    // Define beam direction
    Three_Mom p_hat(0, 0, 1);
    //Three_Mom p_hat=ttbar.Vect().Unit(); // beam axis should be parallel to z!
    
    double y_p = p_hat.Dot(k_hat);

    return abs(y_p);
}
        
double CalculatePoptThetaW(Four_Mom_ref _had_W, Four_Mom_ref _had_top, Four_Mom_ref ttbar, Four_Mom_ref jet1, Four_Mom_ref jet2){
    Four_Mom had_top=_had_top;
    had_top.Boost(-ttbar.BoostVector()); 
    
    Four_Mom had_W=_had_W;
    had_W.Boost(-ttbar.BoostVector());  
    had_W.Boost(-had_top.BoostVector()); 
    
    Four_Mom jet_softer=jet2;// any one is fine. we do abs
    jet_softer.Boost(-ttbar.BoostVector());
    jet_softer.Boost(-had_top.BoostVector());
    jet_softer.Boost(-had_W.BoostVector());
    
    double cos_theta_W=jet_softer.Vect().Unit().Dot(had_W.Vect().Unit());
    return abs(cos_theta_W); 
}

double get_theta_lep(std::string axis, Four_Mom_ref _lep_top, Frame_ref helix, Four_Mom_ref ttbar, Four_Mom_ref par){ 
    Four_Mom lep_top = _lep_top;
    lep_top.Boost(-ttbar.BoostVector());
    
    Four_Mom lep = par;
    lep.Boost(-ttbar.BoostVector());
    lep.Boost(-lep_top.BoostVector());
    
    return helix.at(axis).Dot(lep.Vect().Unit());
}
    
double get_theta_jet_opt(std::string axis, Four_Mom_ref _had_top, Frame_ref helix, Four_Mom_ref ttbar, double cos_theta_W, Four_Mom_ref jet1, Four_Mom_ref jet2){ 
    cos_theta_W=abs(cos_theta_W); // for safe
    Four_Mom had_top = _had_top;
    had_top.Boost(-ttbar.BoostVector());

    Four_Mom jet_harder=jet1;
    Four_Mom jet_softer=jet2;
    jet_softer.Boost(-ttbar.BoostVector());
    jet_softer.Boost(-had_top.BoostVector());
    jet_harder.Boost(-ttbar.BoostVector());
    jet_harder.Boost(-had_top.BoostVector());
    
    if(jet_harder.E()<jet_softer.E()){ // need to judge again!
        Four_Mom jet_temp=jet_softer;
        jet_softer=jet_harder;
        jet_harder=jet_temp;
    }
    Three_Mom q_softer(jet_softer.Vect().Unit());
    Three_Mom q_harder(jet_harder.Vect().Unit());
    
    const double mt=172.76; // GeV
    const double mw=80.377;
    double f=(3./4)*(mt*mt/(mt*mt+2*mw*mw))*(1-cos_theta_W*cos_theta_W) + (3./8)*(2*mw*mw/(mt*mt+2*mw*mw))*(1-abs(cos_theta_W))*(1-abs(cos_theta_W));
    double fn=(3./4)*(mt*mt/(mt*mt+2*mw*mw))*(1-cos_theta_W*cos_theta_W) + (3./8)*(2*mw*mw/(mt*mt+2*mw*mw))*(1+abs(cos_theta_W))*(1+abs(cos_theta_W));
    
    // use simplified ratio as the author used
    //double f=(3./4)*(7./10)*(1-cos_theta_W*cos_theta_W) + (3./8)*(3./10)*(1-abs(cos_theta_W))*(1-abs(cos_theta_W));
    //double fn=(3./4)*(7./10)*(1-cos_theta_W*cos_theta_W) + (3./8)*(3./10)*(1+abs(cos_theta_W))*(1+abs(cos_theta_W));

    Three_Mom p3_opt=(fn/(f+fn))*q_softer + (f/(f+fn))*q_harder;

    return helix.at(axis).Dot(p3_opt.Unit());
}

