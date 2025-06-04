using HCubature  # loaded by CalculusWithJulia
using Dates

rmin=0.0
rmax=6.0


# Parameters
c       =   0.000161
a       =   0.000000
alfa1   =   2.348060
alfa2   =   3.303136
alfa3   =   0.556262
alfa12  =   0.209718
c12     =   0.000000
d12     =   0.000000
alfa13  =   0.037843
alfa23  =   0.011527    

       

global pars=[c,a,alfa1,alfa2,alfa3,alfa12,c12,d12,alfa13,alfa23]
global Z=2

function psi(rho1,phi1,z1,rho2,phi2,z2,rho3,phi3,z3)
      a    = pars[2]
      al1  = pars[3]
      al2  = pars[4]
      al3  = pars[5]
      al12 = pars[6]
      c12  = pars[7]
      d12  = pars[8]
      al13 = pars[9]
      al23 = pars[10]
      r1=(rho1^2+z1^2)
      r2=(rho2^2+z2^2)
      r3=(rho3^2+z3^2)
      r12 = (rho1^2+rho2^2−2*rho1*rho2*cos(phi1−phi2)+(z1−z2)^2)^(1/2)
      r13 = (rho1^2+rho3^2−2*rho1*rho3*cos(phi1−phi3)+(z1−z3)^2)^(1/2)
      r23 = (rho2^2+rho3^2−2*rho2*rho3*cos(phi2−phi3)+(z2−z3)^2)^(1/2)
      wf=z3*(1+a*r12)*exp(-al1*r1-al2*r2-al3*r3+al12*r12*(1+c12*r12)/(1+d12*r12)+al13*r13+al23*r23)
      return wf
      end


function psid(rho1,phi1,z1,rho2,phi2,z2,rho3,phi3,z3,iflag)
      a    = pars[2]
      al1  = pars[3]
      al2  = pars[4]
      al3  = pars[5]
      al12 = pars[6]
      c12  = pars[7]
      d12  = pars[8]
      al13 = pars[9]
      al23 = pars[10]
      r1=(rho1^2+z1^2)
      r2=(rho2^2+z2^2)
      r3=(rho3^2+z3^2)
      r12 = (rho1^2+rho2^2−2*rho1*rho2*cos(phi1−phi2)+(z1−z2)^2)^(1/2)
      r13 = (rho1^2+rho3^2−2*rho1*rho3*cos(phi1−phi3)+(z1−z3)^2)^(1/2)
      r23 = (rho2^2+rho3^2−2*rho2*rho3*cos(phi2−phi3)+(z2−z3)^2)^(1/2)
      expfac=exp(-al1*r1-al2*r2-al3*r3+al12*r12*(1+c12*r12)/(1+d12*r12)+al13*r13+al23*r23)
      pfac = z3*(1+a*r12)
      if (iflag==1)
        wfm=al1*pfac*expfac
      end
      if (iflag==2)
        wfm=al2*pfac*expfac
      end
      if (iflag==3)
        wfm=(al3*pfac+0.0*(1+a*r12))*expfac
      end
      return wfm
      end

function denominator(rho2,phi2,z2,rho3,phi3,z3)
         rho1=0.0
	 phi1=0.0
	 z1=0.0
	 c=pars[1]
	 psi123=psi(rho1,phi1,z1,rho2,phi2,z2,rho3,phi3,z3)
	 psi213=psi(rho2,phi2,z2,rho1,phi1,z1,rho3,phi3,z3)
	 psi321=psi(rho3,phi3,z3,rho2,phi2,z2,rho1,phi1,z1)
	 psi132=psi(rho1,phi1,z1,rho3,phi3,z3,rho2,phi2,z2)
	 psi312=psi(rho3,phi3,z3,rho1,phi1,z1,rho2,phi2,z2)
	 psi231=psi(rho2,phi2,z2,rho3,phi3,z3,rho1,phi1,z1)
	 P1 = (1/(2*3^(1/2)))*(2*psi123+2*psi213-psi321-psi132-psi231-psi312)
	 P2 = (1/2)*(psi321-psi132-psi231+psi312)
	 P3 = (1/2)*(psi321-psi132+psi231-psi312)
	 P4 = (1/(2*3^(1/2)))*(2*psi123-2*psi213+psi321+psi132-psi231-psi312)
	 psit2=(P1 + c*P3)^2 + (P2 + c*P4)^2
	 jac=rho2*rho3
	 res = psit2*jac
	 return res 
end

function numerator(rho2,phi2,z2,rho3,phi3,z3)
         rho1=0.0
	 phi1=0.0
	 z1=0.0
	 c=pars[1]
	 psi123=psi(rho1,phi1,z1,rho2,phi2,z2,rho3,phi3,z3)
	 psi213=psi(rho2,phi2,z2,rho1,phi1,z1,rho3,phi3,z3)
	 psi321=psi(rho3,phi3,z3,rho2,phi2,z2,rho1,phi1,z1)
	 psi132=psi(rho1,phi1,z1,rho3,phi3,z3,rho2,phi2,z2)
	 psi312=psi(rho3,phi3,z3,rho1,phi1,z1,rho2,phi2,z2)
	 psi231=psi(rho2,phi2,z2,rho3,phi3,z3,rho1,phi1,z1)
	 #----------------------------
	 psi123d=psid(rho1,phi1,z1,rho2,phi2,z2,rho3,phi3,z3,1)
	 psi213d=psid(rho2,phi2,z2,rho1,phi1,z1,rho3,phi3,z3,2)
	 psi321d=psid(rho3,phi3,z3,rho2,phi2,z2,rho1,phi1,z1,3)
	 psi132d=psid(rho1,phi1,z1,rho3,phi3,z3,rho2,phi2,z2,1)
	 psi312d=psid(rho3,phi3,z3,rho1,phi1,z1,rho2,phi2,z2,2)
	 psi231d=psid(rho2,phi2,z2,rho3,phi3,z3,rho1,phi1,z1,3)
	 #
	 P1 = (1/(2*3^(1/2)))*(2*psi123+2*psi213-psi321-psi132-psi231-psi312)
	 P2 = (1/2)*(psi321-psi132-psi231+psi312)
	 P3 = (1/2)*(psi321-psi132+psi231-psi312)
	 P4 = (1/(2*3^(1/2)))*(2*psi123-2*psi213+psi321+psi132-psi231-psi312)
	 P1d = (1/(2*3^(1/2)))*(2*psi123d+2*psi213d-psi321d-psi132d-psi231d-psi312d)
	 P2d = (1/2)*(psi321d-psi132d-psi231d+psi312d)
	 P3d = (1/2)*(psi321d-psi132d+psi231d-psi312d)
	 P4d = (1/(2*3^(1/2)))*(2*psi123d-2*psi213d+psi321d+psi132d-psi231d-psi312d)
	 #----------------------------------------------------------
	 psit2=(P1 + c*P3)*(P1d + c*P3d) + (P2 + c*P4)*(P2d + c*P4d)
	 jac=rho2*rho3
	 res = psit2*jac
	 return res 
end




denominator(v) = denominator(v...)  # denominator accepts a vector
numerator(v) = numerator(v...)  # denominator accepts a vector


println("Expectation values of the Lithium-like atoms 1s^22p state")
println("version updated:", now())

eps=0.1
Nmax = 10000000
den=hcubature(denominator,(rmin,0,-rmax,rmin,0,-rmax),(rmax,2*pi,rmax,rmax,2*pi,rmax),maxevals=Nmax)
println("denominator                :",den)
num=hcubature(numerator,(rmin,0,-rmax,rmin,0,-rmax),(rmax,2*pi,rmax,rmax,2*pi,rmax),maxevals=Nmax)
println("numerator                :",num)
println("Cusp_ee                  :",num[1]/den[1])


      
      
   
          
       
                     
       
        
       
            
            
            
            
            
            
            
            
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    

