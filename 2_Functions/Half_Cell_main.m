  function J = Half_Cell_main(x,data) 
  battery_type=data.battery_type;
  % Due to hysterisis, Un for charge and discharge are different
  % Those are four stochimotery parameters for discharge
        thetap100=x(1);
        thetap0=x(2);
        thetan100=x(3);
        thetan0=x(4);
        
        if battery_type==1
  % Those are four stochimotery parameters for charge      
        thetap100ch=x(5);
        thetap0ch=x(6);
        thetan100ch=x(7);
        thetan0ch=x(8);
        end

        warning off
        try

        Vocfull=data.Vocfull_DC;
        
        totnum=data.intnum_DC;
        intsocp=linspace(thetap100,thetap0,totnum);
        intOCVp=interp1(data.OCV.pe_socdis,data.OCV.pe_Updis,intsocp,'linear','extrap');
        
        intsocn=linspace(thetan0,thetan100,totnum);
        intOCVn=interp1(data.OCV.ne_socdis,data.OCV.ne_Undis,intsocn,'linear','extrap');
        
        Vsim=intOCVp-flip(intOCVn);
        Vexp=Vocfull;
  
        J1=rms(Vsim-Vexp)*1000;      
        J=J1;
        if battery_type==1
        Vocfull=data.Vocfull_CC;
        
        totnum=data.intnum_CC;
        intsocp=linspace(thetap100ch,thetap0ch,totnum);
        intOCVp=interp1(data.OCV.pe_socch,data.OCV.pe_Upch,intsocp,'linear','extrap');
        
        intsocn=linspace(thetan0ch,thetan100ch,totnum);
        intOCVn=interp1(data.OCV.ne_socch,data.OCV.ne_Unch,intsocn,'linear','extrap'); 
        
        Vsim=intOCVp-flip(intOCVn);
        Vexp=Vocfull;
            J2=rms(Vsim-Vexp)*1000;           
            J=J1+J2;
        end
        catch
            J = 10^3;
        end
            
    end