function NACA4points(number)
        #TODO: This function comes from the ME 415 repository
        naca = "$number"
        n = 80
        # parse string
        e = parse(Int64, naca[1])/100.0
        p = parse(Int64, naca[2])/10.0
        t = parse(Int64, naca[3:4])/100.0
    
        # cosing spacing
        theta = range(0, pi, length=n)
        x = (1.0 .- cos.(theta))/2.0
    
        T = @. 10*t*(0.2969*sqrt(x) - 0.126*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)
    
        ybar = similar(x)
        for i = 1:n
            if x[i] < p # was <=, but causes problems at the leading edge of symmetric airfoils
                ybar[i] = e/p^2 * (2*p*x[i] - x[i]^2)
            else
                ybar[i] = e/(1-p)^2 * (1 - 2*p + 2*p*x[i] - x[i]^2)
            end
        end
    
        yu = ybar + T/2.0
        yl = ybar - T/2.0
    
        y = [yu[end:-1:1]; yl]
        x = [x[end:-1:1]; x]
    
    return x,y
end

function NACA6points(number)
    error("NACA6 not implimented")
    return x,y
end