module fem

export stiffMat
export loadVct

    function stiffMat(ea, ei, el, n, v)
        # streching part
        ma = [
            1/el -1/el;
            -1/el 1/el
        ]
        mn = (n + 1)*3
        ukn = zeros(mn, mn)

        for i = 1:n
            alpha = v[i]
            ks = ma*ea
            # bending part
            mi = [
                12/el^3 6/el^2 -12/el^3 6/el^2;
                6/el^2 4/el -6/el^2 2/el;
                -12/el^3 -6/el^2 12/el^3 -6/el^2;
                6/el^2 2/el -6/el^2 4/el
            ]
            kb = alpha*mi*ei

            uk0 = zeros(6, 6)
            uk0[1:3:4, 1:3:4] = ks
            uk0[2:3, 2:3] = kb[1:2, 1:2]
            uk0[2:3, 5:6] = kb[1:2, 3:4]
            uk0[5:6, 2:3] = kb[3:4, 1:2]
            uk0[5:6, 5:6] = kb[3:4, 3:4]
            pa = 1 + 3*(i-1)
            pb = 6 + 3*(i-1)
            ukn[pa:pb, pa:pb] = ukn[pa:pb, pa:pb] + uk0
        end

        # boundary conditions
        for i in [0, n]
            p = [ 1+3*i 2+3*i 3+3*i ]
            ukn[p[1]:p[3], :] = zeros(3, mn)
            ukn[:, p[1]:p[3]] = zeros(mn, 3)
            for j in p
                ukn[j, j] = 1
            end
        end
        return K = ukn
    end

    function loadVct(q, el, n)
        # load
        vl = [
                0; q*el/2; q*el^2/12;
                0; q*el/2; q*el^2/12
        ]
        vn = (n + 1)*3
        f0 = zeros(vn, 1)
        for i = 1:n
                pa = 1 + 3*(i-1)
                pb = 6 + 3*(i-1)
                f0[pa:pb, :] = f0[pa:pb, :] + vl
        end
        # boundary conditions
        for i in [0, n]
            p = [ 1+3*i 2+3*i 3+3*i ]
            f0[p[1]:p[3], :] = zeros(3, 1)
        end
        return f = f0
    end
end  # module fem
