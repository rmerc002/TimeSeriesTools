function [CP, MPposneg, MPpospos, Tposnew] = ContrastProfileI_old(Tneg,Tpos, t, MPposneg, MPpospos, m)
    Tposnew = [Tpos, t];
    Spos = Tposnew(end-m+1:end);
    
    MPposneg(end+1) = min(MASS_V2(Tneg, Spos));
    
    DPpospos = MASS_V2(Tpos, Spos);
    MPpospos = min(MPpospos, DPpospos);
    MPpospos(end+1) = min(DPpospos);
    
    CP = MPposneg - MPpospos;
end