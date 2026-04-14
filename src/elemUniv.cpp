#include "elemUniv.h"

void ElemUniv::funkcjaKsztaltu(unsigned int npc) 
{
    unsigned int totalPC = npc * npc;
    GaussLegendreData gauss(npc);
    unsigned int index = 0;

    for (unsigned int i = 0; i < npc; i++) 
    {
        for (unsigned int j = 0; j < npc; j++) 
        {
            double ksi, eta;
            if (npc == 2)
            {
                ksi = gauss.nodes[j];
                eta = gauss.nodes[i];
            }
            else
            {
                ksi = gauss.nodes[i];
                eta = gauss.nodes[j];
            }

            // ksi = gauss.nodes[i];
            // eta = gauss.nodes[j];

            dNdKsi[index][0] = -0.25 * (1 - eta);
            dNdKsi[index][1] =  0.25 * (1 - eta);
            dNdKsi[index][2] =  0.25 * (1 + eta);
            dNdKsi[index][3] = -0.25 * (1 + eta);

            dNdEta[index][0] = -0.25 * (1 - ksi);
            dNdEta[index][1] = -0.25 * (1 + ksi);
            dNdEta[index][2] =  0.25 * (1 + ksi);
            dNdEta[index][3] =  0.25 * (1 - ksi);

            N[index][0] = 0.25 * (1 - ksi) * (1 - eta);
            N[index][1] = 0.25 * (1 + ksi) * (1 - eta);
            N[index][2] = 0.25 * (1 + ksi) * (1 + eta);
            N[index][3] = 0.25 * (1 - ksi) * (1 + eta);

            index++;
        }
    }
    for (unsigned int gp = 0; gp < npc; ++gp)
    {
        double ksi, eta;

        ksi = gauss.nodes[gp]; 
        eta = -1.0;
        s[0].N[gp][0] = 0.25 * (1 - ksi) * (1 - eta);
        s[0].N[gp][1] = 0.25 * (1 + ksi) * (1 - eta);
        s[0].N[gp][2] = 0.25 * (1 + ksi) * (1 + eta);
        s[0].N[gp][3] = 0.25 * (1 - ksi) * (1 + eta);

        ksi = 1.0;
        eta = gauss.nodes[gp];
        s[1].N[gp][0] = 0.25 * (1 - ksi) * (1 - eta);
        s[1].N[gp][1] = 0.25 * (1 + ksi) * (1 - eta);
        s[1].N[gp][2] = 0.25 * (1 + ksi) * (1 + eta);
        s[1].N[gp][3] = 0.25 * (1 - ksi) * (1 + eta);

        ksi = gauss.nodes[gp];
        eta = 1.0;
        s[2].N[gp][0] = 0.25 * (1 - ksi) * (1 - eta);
        s[2].N[gp][1] = 0.25 * (1 + ksi) * (1 - eta);
        s[2].N[gp][2] = 0.25 * (1 + ksi) * (1 + eta);
        s[2].N[gp][3] = 0.25 * (1 - ksi) * (1 + eta);

        ksi = -1.0;
        eta = gauss.nodes[gp];
        s[3].N[gp][0] = 0.25 * (1 - ksi) * (1 - eta);
        s[3].N[gp][1] = 0.25 * (1 + ksi) * (1 - eta);
        s[3].N[gp][2] = 0.25 * (1 + ksi) * (1 + eta);
        s[3].N[gp][3] = 0.25 * (1 - ksi) * (1 + eta);
    }
}

ElemUniv::ElemUniv(unsigned int npc)
{
    unsigned int totalPC = npc * npc;
    dNdKsi = new double*[totalPC];  
    dNdEta = new double*[totalPC];
    N = new double*[totalPC];
    for (unsigned int i = 0; i < totalPC; ++i) 
    {
        dNdKsi[i] = new double[4];
        dNdEta[i] = new double[4];
        N[i] = new double[4];
    }
    for (int i = 0; i < 4; ++i) 
    {
        s[i].N = new double*[npc];
        for (unsigned int j = 0; j < npc; ++j) 
        {
            s[i].N[j] = new double[4];
        }
    }
    funkcjaKsztaltu(npc);
}