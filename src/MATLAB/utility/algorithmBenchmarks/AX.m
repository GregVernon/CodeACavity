function ax = AX(AX_jm1i,AX_jp1i,AX_jim1,AX_jip1,AX_ji)
    ax = AX_jm1i + AX_jp1i + AX_jim1 + AX_jip1 - 4*AX_ji;
end