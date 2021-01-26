#We  are using this to make the Cd v. Velocity Plots
#import necessary libraries

using Plots


function air_Density_Calculation(h)
    temp = 59 - .00356 * h;         #h = altitude(ft), temp = Fereinheidt.

    rho = 2116 * [(temp + 459.7)/ 518.6]^5.256;
    #Density as a function of tempurature. Only of temp below 36152 ft.

    return rho;

end
function root_Mean_Chord(c_r, c_t)
    c_mac = (2/3) * (c_r + c_t + ((c_r*c_t)/(c_r+c_t)));
    return c_mac
end

function form_Factor_Calc(t_c, c_mac, Delta)
    #requires tip chord, root mean Chord, and sweep angle (Delta)
    #returns value of k

    Z = 2*cos(Delta)
    k = 1 + Z*(t_c/c_mac) + 100*(t_c/c_mac)^4

end
#Function Needs: Speed, Altitude, Surface Area, Tip Chord, Root Mean Chord, Sweep Angle, and mu.
function parasitic_Drag_Calculation(V_inf, q_inf, S_wet, tip2chord, c_mac, Delta, mu)
    #S_Wet = 2*(1+ 0.2*(tip2chord)) * S_wing;  #S_Wing passed in vector.
    #q_inf = (1/2) * p_inf * (V_inf^2);
    k = form_Factor_Calc(t, c, Delta);
    #Pass root chord and tip chord lengths, Only needed for Cessna, so don't do it here.
    RE = (p_inf * V_inf * c)/(mu);
    C_f = 1.328/sqrt(RE);           #Laminar flow, we won't worry about turbulent flow.

    #Use atmospheric caluculator online to find specific values for mu, T, etc.

    #mu is passed on vector, c is the mean aerodynamic chord of the aircraft.
    #mu is the dynamic viscocity of the fluid.

    #fineness ratio, fr, is the ratio of an airplanes length to its width.c

    #D_p = k * C_f * q_inf * S_Wet
    #q_inf = (1/2) * p_inf * (V_inf)^2
    #RE = (p*V*c) / (mu)
    #   p = fluid density
    #   mu = fluid viscocity

    D_p = k * C_f * q_inf * S_Wet;  #So far only need k and C_f

    return D_p;
end
function induced_Drag_Calculation(L, v_inf, h, )
    p_inf

    D_i = L^2 / (q_inf * pi * b^2 * e_inv)
end
function main()
    #We need 2 values, induced drag and parasitic drag.
    #If plane cruises above Mach 0.3 (>230 mph), include Compressiblity Drag

    #Choose 2 planes whose values you mimic. One should be glider esq. (high aspect ratio) and the other more like a fighter jet (low aspect ratio)

    #Cessna 172R
        # Aspect Ratio: 7.32 (according to wikipedia)
        # Length: 27 ft, 2 in
        # Wing Span: 36 ft, 1 in.
        # Surface Area: 174 ft^2, S_wing
        # Aspect Ratio: 4.24153... (according to my calculations)
        # Cruising Speed: 140 mph (top 188 mph)
        # Stalling speed: 54 mph
        # Loaded Wight: 2450 lb
        # Cruising Alt: 6000 ft., mu = 3.6162e-7   lb-s / ft^2
        # Overall Length: 26 ft, 11 in.


    #Vought V-173, Flying Pancake (oldmachinepress.com)
        # Span: 23 ft, 4 in.
        # Surface Area: 427 ft^2, S_wing
        # Aspect Ratio = 1.275....
        # Max Speed: 138 mph
        # Cruising Speed: 75 mph
        # Stalling Speed:
        # Cruising Altitude: 28000 ft., mu = 3.1510e-7  lb-s / ft^2
        # Root mean chord: SurfaceArea / Span

        # Loaded Weight = 3050 lb

        #Calculate the root mean chord, c, before calling drag functions.

        rho = air_Density_Calculation(h);
        S_Wet = 2*(1+ 0.2*(tip2chord)) * S_wing;  #S_Wing passed in vector.
        q_inf = (1/2) * p_inf * (V_inf^2);

        mu = 3.1510e-7; #28000 ft.
        # mu = 3.6162e-7; #6000 ft.

    #Parasitic Drag:
        D_p = parasitic_Drag_Calculation(velocity, rho, S_Wet, tip2chord, c_mac, delta, mu);
        #D_p = k * C_f * q_inf * S_wet
        #k = Pressure drag constant
        #C_f = Skin Friction drag constant
        #q_inf = Wind speed over wing
        #S_wet = 2*(1+ 0.2*(t/c)) * S_wing
        # c_mean is the root mean chord.
        # t_c is the thickness to chord ratio, as a decimal-percent

        #Function Needs: Speed, Altitude, Surface Area, Tip Chord, Root Mean Chord, Sweep Angle, and mu.

    #Induced Drag
        df = 1;                         #Diameter of the fuselage.
        e_inv = 0.98 *(1-(2*(d_f/b)^2))
        D_i = induced_Drag_Calcuation();
        #D_i = L^2 / (q_inf * pi * b^2 * e_inv)
        #C_di = C_l^2 / (pi * AR * e_inv)

        #e_inv(Invicide span efficiency) = 0.98 * [1-2*(d_f/b)^2]
        #d_f = diameter of the fuselage
        #L is lift, which is just equal to the weightof the aircraft.

        #D = D_p + D_i

        #Velocity range (0, cruising speed), Do this to start
        #Find the point at which the two drags equalize, this is the minimum
        # value we want
end
