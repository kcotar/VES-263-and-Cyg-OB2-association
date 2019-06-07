import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, join, unique, hstack, vstack
from os import path
from scipy.stats import multivariate_normal, norm
import astropy.coordinates as coord
import astropy.units as un
import imp
import joblib
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN, KMeans
from gaia_data_queries import get_data_subset

plt.rcParams['font.size'] = 15
# plt.rcParams['font.family'] = 'cursive'
# plt.rcParams['font.cursive'] = 'Textile'

d1 = Table.read('Berlanas_2018.csv', format='ascii.csv')
d1 = d1[d1['Sample'] == 'Cygnus OB2']
d2 = Table.read('Comeron_2012.csv', format='ascii.csv')
d2 = d2[d2['Field'] == 'Cygnus OB2']
d3 = Table.read('Kiminki_2007.csv', format='ascii.csv')
# d4 = Table.read('Kharchenko_2004.csv', format='ascii.csv')
# d4 = d4[d4['Seq'] == 490]
# d4.write('Kharchenko_2004_Cyg_OB2.csv', format='ascii.csv')
d4 = Table.read('Kharchenko_2004_Cyg_OB2.csv', format='ascii.csv')
d5 = Table.read('cyg_ob2_simbad.csv', format='ascii.csv')
d6 = Table.read('Wright_2015.csv', format='ascii.csv')

ves_263_source = 2067735760399789568

cyg_ob2 = vstack((
                  d1['source_id', 'ra', 'dec'],
                  d2['source_id', 'ra', 'dec'],
                  # d3['source_id', 'ra', 'dec'],
                  # d4['source_id', 'ra', 'dec'],
                  # d5['source_id', 'ra', 'dec'],
                  d6['source_id', 'ra', 'dec']
                  ))

#print cyg_ob2['source_id', 'ra', 'dec', 'parallax', 'pmra', 'pmdec']
cyg_ob2 = unique(cyg_ob2, keys='source_id')

print 'Members:', len(cyg_ob2), len(np.unique(cyg_ob2['source_id']))
print 'Number VES 263:', np.sum(cyg_ob2['source_id'] == ves_263_source)
print 'Number VES 263:', np.sum(d1['source_id'] == ves_263_source), np.sum(d2['source_id'] == ves_263_source), np.sum(d3['source_id'] == ves_263_source), np.sum(d4['source_id'] == ves_263_source), np.sum(d5['source_id'] == ves_263_source), np.sum(d6['source_id'] == ves_263_source)
print 'Coords median:', np.nanmedian(cyg_ob2['ra']), np.nanmedian(cyg_ob2['dec'])

gaia_file = 'cyg_ob2_gaia.csv'
if path.isfile(gaia_file):
    gaia_data_all = Table.read(gaia_file, format='ascii.csv')
else:
    print 'Gaia DR2 query sent'
    gaia_data_all = get_data_subset(np.nanmedian(cyg_ob2['ra']),
                                np.nanmedian(cyg_ob2['dec']),
                                3., None, dist_span=None, rv_only=False)
    gaia_data_all.write(gaia_file, format='ascii.csv')

for col in ['pmra', 'pmdec', 'parallax', 'rv']:
    gaia_data_all['rv'][gaia_data_all['rv'].filled(-1e6) == -1e6] = np.nan

print 'Gaia DR2 in cone:', len(gaia_data_all)
cyg_ob2 = gaia_data_all[np.in1d(gaia_data_all['source_id'], cyg_ob2['source_id'])]

gaia_ra_dec = coord.ICRS(ra=gaia_data_all['ra']*un.deg, dec=gaia_data_all['dec']*un.deg)
ob2_ra_dec = coord.ICRS(ra=np.nanmedian(cyg_ob2['ra'])*un.deg, dec=np.nanmedian(cyg_ob2['dec'])*un.deg)
idx_use_coord = gaia_ra_dec.separation(ob2_ra_dec) <= 1.*un.deg

plt.figure(figsize=(7., 3.5))
plt.hist(cyg_ob2['parallax'], range=(0.4, 0.85),
         bins=np.arange(0.36, 0.92, 0.02),
         color='black', histtype='step', lw=1.5)
plt.axvline(0.5962, c='red')
plt.axvline(0.58, c='black', ls='--', alpha=0.8)
plt.axvline(0.68, c='black', ls='--', alpha=0.5)
plt.axvline(0.48, c='black', ls='--', alpha=0.5)
plt.xlim(0.39, 0.86)
plt.ylabel('Number of stars in a bin')
plt.xlabel('Parallax (mas)')
plt.gca().tick_params(axis="y", direction="in")
plt.gca().tick_params(axis="x", direction="in")
# plt.grid(ls='--', c='black', alpha=0.2)
plt.tight_layout()
plt.savefig('cygob2_parallax.pdf', dpi=250)
plt.close()

# # get all the data and perform clusterings
# ob2_data = gaia_data_all[np.in1d(gaia_data_all['source_id'], cyg_ob2['source_id'])]
# ves263_data = gaia_data_all[np.in1d(gaia_data_all['source_id'], ves_263_source)]
# ob2_data_cluster = ob2_data['pmra', 'pmdec', 'parallax'].to_pandas().values
# # standardizatrion
# ob2_data_cluster = (ob2_data_cluster - np.median(ob2_data_cluster, axis=0)) / np.std(ob2_data_cluster, axis=0)
# # dbscan
# # labels = DBSCAN(0.5, min_samples=5).fit_predict(ob2_data_cluster)
# labels = KMeans(8, n_init=10, n_jobs=50).fit_predict(ob2_data_cluster)
# plt.scatter(ob2_data['ra'], ob2_data['dec'], c=labels, lw=0, s=15, cmap='Set1')
# plt.colorbar()
# plt.scatter(ves263_data['ra'], ves263_data['dec'], c='black', lw=0, s=20, marker='X')
# plt.show()
# plt.close()
#
# tsne_2d = TSNE(n_components=2, perplexity=10.0, method='barnes_hut', angle=0.5).fit_transform(ob2_data_cluster)
# plt.scatter(tsne_2d[:, 0], tsne_2d[:, 1], c=labels, lw=0, s=10, cmap='Set1')
# plt.show()
# plt.close()
#
# for col in ['ra', 'dec', 'pmra', 'pmdec', 'parallax']:
#     plt.scatter(tsne_2d[:, 0], tsne_2d[:, 1], c=ob2_data[col], lw=0, s=10)
#     plt.title(col)
#     plt.show()
#     plt.close()
#
# raise SystemExit

d_parallax = 0.1
# mean_parallax = 0.3
for mean_parallax in [0.58]:  # [0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15]:
    idx_use = np.logical_and(idx_use_coord, np.abs(gaia_data_all['parallax'] - mean_parallax) <= d_parallax)

    gaia_data = gaia_data_all[idx_use]

    gaia_ra_dec = coord.ICRS(ra=gaia_data['ra']*un.deg, dec=gaia_data['dec']*un.deg, distance=1e3/gaia_data['parallax']*un.pc)
    gaia_gal = gaia_ra_dec.transform_to(coord.Galactocentric())

    idx_ob2 = np.in1d(gaia_data['source_id'], cyg_ob2['source_id'])
    idx_ves263 = gaia_data['source_id'] == ves_263_source

    idx_ob2 = np.logical_xor(idx_ob2, idx_ves263)

    ves_263 = gaia_data[gaia_data['source_id'] == ves_263_source]
    print 'VES:', ves_263['ra'], ves_263['dec'], ves_263['parallax'],  ves_263['parallax_error']

    # plt.scatter(gaia_gal.x, gaia_gal.y, lw=0, s=6, alpha=0.2, c='black')
    # plt.scatter(gaia_gal.x[idx_ob2], gaia_gal.y[idx_ob2], lw=0, s=12, alpha=1., c='red')
    # plt.scatter(gaia_gal.x[idx_ves263], gaia_gal.y[idx_ves263], lw=0, s=50, alpha=1., c='blue')
    # plt.show()
    # plt.close()

    idx_inparams = np.logical_and(gaia_data['pmra'] >= -0.5, gaia_data['pmra'] <= 0.5)
    idx_inparams = np.logical_and(idx_inparams, gaia_data['pmdec'] >= -4)
    idx_inparams = np.logical_and(idx_inparams, gaia_data['pmdec'] <= -2.5)
    idx_inparams = np.logical_and(idx_inparams, gaia_data['parallax'] >= 0.65)
    idx_inparams = np.logical_and(idx_inparams, gaia_data['parallax'] <= 0.85)
    print 'In selection:', np.sum(idx_inparams)
    print 'All stars:', len(gaia_data)
    print 'OB2 stars:', np.sum(idx_ob2)

    # plt.scatter(gaia_data['ra'], gaia_data['dec'], lw=0, s=6, alpha=0.2, c='black')
    # plt.scatter(gaia_data['ra'][idx_inparams], gaia_data['dec'][idx_inparams], lw=0, s=6, alpha=1., c='green')
    # plt.scatter(gaia_data['ra'][idx_ob2], gaia_data['dec'][idx_ob2], lw=0, s=12, alpha=1., c='red')
    # plt.scatter(ves_263['ra'], ves_263['dec'], lw=0, s=50, alpha=1., c='blue')
    # plt.quiver(gaia_data['ra'][idx_ob2], gaia_data['dec'][idx_ob2],
    #            gaia_data['pmra'][idx_ob2] - np.nanmedian(gaia_data['pmra']),
    #            gaia_data['pmdec'][idx_ob2] - np.nanmedian(gaia_data['pmdec']),
    #            color='green', width=0.003)
    # plt.xlabel('RA')
    # plt.ylabel('DEC')
    # plt.show()
    # plt.close()
    #
    # plt.scatter(gaia_data['phot_bp_mean_mag'][idx_inparams]-gaia_data['phot_rp_mean_mag'][idx_inparams],
    #             gaia_data['phot_g_mean_mag'][idx_inparams] - 2.5*np.log10(((1e3/gaia_data['parallax'][idx_inparams])/10.)**2)
    #             )
    # plt.gca().invert_yaxis()
    # plt.show()
    # plt.close()

    # pmra/pmdec density estimation plot
    pmra_vals = np.linspace(-12, 9, 350)
    pmdec_vals = np.linspace(-14, 4, 350)
    xx, yy = np.meshgrid(pmra_vals, pmdec_vals)
    pos = np.empty(xx.shape + (2,))
    pos[:, :, 0] = xx
    pos[:, :, 1] = yy

    file_density = 'density_{:.2f}_{:.2f}.pkl'.format(mean_parallax, d_parallax)
    if path.isfile(file_density):
        total_density = joblib.load(file_density)
    else:
        total_density = np.full_like(xx, fill_value=0.)
        for i_s, star in enumerate(gaia_data):
            if i_s % 500 == 0:
                print i_s
            cov = np.array([[star['pmra_error'], 0], [0, star['pmdec_error']]])
            mean = np.array([star['pmra'], star['pmdec']])
            total_density += multivariate_normal.pdf(pos, mean=mean, cov=cov)
        joblib.dump(total_density, file_density)

    plt.imshow(total_density, interpolation='none', origin='lower', cmap='viridis',
               extent=[-12, 9, -14, 4],
               vmin=0,  # np.percentile(total_density, 0),
               vmax=1300,  # np.percentile(total_density, 100)
               )
    plt.colorbar()
    plt.contour(xx, yy, total_density, [0, 50, 200, 350, 500, 650, 800, 950, 1100, 1250, 1400, 1550, 1700],
                colors='black', linewidths=0.5)
    # plt.scatter(gaia_data['pmra'], gaia_data['pmdec'], lw=0, s=1, alpha=0.2, c='black')
    # plt.scatter(gaia_data['pmra'][idx_ob2], gaia_data['pmdec'][idx_ob2], lw=0, s=1, alpha=1., c='C0')
    plt.scatter(ves_263['pmra'], ves_263['pmdec'], lw=0, s=2, alpha=1., c='red')
    plt.title('Number of stars: {:.0f}'.format(len(gaia_data)))
    plt.xlim(np.min(pmra_vals), np.max(pmra_vals))
    plt.ylim(np.min(pmdec_vals), np.max(pmdec_vals))
    plt.gca().tick_params(axis="y", direction="in")
    plt.gca().tick_params(axis="x", direction="in")
    # plt.show()
    plt.tight_layout()
    plt.savefig('pmra_pmdec_dist_{:.2f}_{:.2f}.pdf'.format(mean_parallax, d_parallax), dpi=350)
    plt.close()

    plt.figure(figsize=(7., 6.))
    plt.scatter(gaia_data['pmra'], gaia_data['pmdec'], lw=0, s=3., alpha=0.5, c='0.6',
                label='', zorder=1)
    # plt.title('Number of stars: {:.0f}'.format(len(gaia_data)))
    # plt.contour(xx, yy, total_density,
    #             [0, 50, 200, 350, 500, 650, 800, 950, 1100, 1250, 1400, 1550, 1700],
    #             colors='black', linewidths=0.45, alpha=0.7, zorder=2)
    plt.contour(xx, yy, total_density,
                [0, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700],
                colors='black', linewidths=0.4, alpha=0.5, zorder=2)
    plt.scatter(gaia_data['pmra'][idx_ob2], gaia_data['pmdec'][idx_ob2], lw=0, s=13, alpha=1., c='green',
                label='Cyg OB2', zorder=3)
    plt.scatter(ves_263['pmra'], ves_263['pmdec'], lw=0, s=70, alpha=1., c='red',
                label='VES 263', marker='X', zorder=3)
    plt.xlim(-6, 2)
    plt.ylim(-8, -1)
    plt.xticks([-5, -3, -1, 1], ['-5', '-3', '-1', '1'])
    plt.yticks([-2, -4, -6, -8], ['-2', '-4', '-6', '-8'])
    plt.gca().tick_params(axis="y", direction="in")
    plt.gca().tick_params(axis="x", direction="in")
    plt.xlabel(u'pmRA (mas yr$^{-1}$)')
    plt.ylabel(u'pmDec (mas yr$^{-1}$)')
    plt.legend(loc=4, framealpha=0.7)
    plt.tight_layout()
    plt.savefig('pmra_pmdec_points_{:.2f}_{:.2f}.pdf'.format(mean_parallax, d_parallax), dpi=350)
    plt.close()

    if mean_parallax == 0.58:
        fig, ax = plt.subplots(2, 1, figsize=(7., 6.))

        ax[0].errorbar(gaia_data['ra'][idx_ob2], gaia_data['pmra'][idx_ob2], yerr=gaia_data['pmra_error'][idx_ob2],
                       elinewidth=0.8, marker='.', ms=7, mew=0, alpha=0.8, c='green', ls='None', label='Cyg OB2',
                       capsize=0, capthick=3)
        ax[0].errorbar(ves_263['ra'], ves_263['pmra'], yerr=ves_263['pmra_error'],
                       elinewidth=0.8, marker='X', ms=8, mew=0, alpha=0.8, c='red', ls='None', label='VES 263',
                       capsize=0, capthick=5)
        lin_fit = np.polyfit(gaia_data['ra'][idx_ob2], gaia_data['pmra'][idx_ob2],
                             deg=1, w=1. / gaia_data['pmra_error'][idx_ob2])
        idx_fit_use = np.logical_and(idx_ob2,
                                     (np.abs(np.polyval(lin_fit, gaia_data['ra']) - gaia_data['pmra'])) <= 1.)
        lin_fit = np.polyfit(gaia_data['ra'][idx_fit_use], gaia_data['pmra'][idx_fit_use],
                             deg=1, w=1. / gaia_data['pmra_error'][idx_fit_use])
        x_lin = np.linspace(306., 310., 25)
        ax[0].plot(x_lin, np.polyval(lin_fit, x_lin), lw=1, ls='--', alpha=0.75, c='black')
        ax[0].set(xlabel=u'RA (deg)', ylabel=u'pmRA (mas yr$^{-1}$)', xlim=(307.1, 308.8), ylim=(-3.8, -1.6),
                  xticks=[307.5, 308, 308.5], xticklabels=['307.5', '308', '308.5'],
                  yticks=[-2, -3], yticklabels=['-2', '-3'])
        ax[0].grid(ls='--', c='black', alpha=0.2)
        ax[0].tick_params(axis="y", direction="in")
        ax[0].tick_params(axis="x", direction="in")
        ax[0].legend(loc=2)

        # plt.scatter(gaia_data['ra'][idx_ob2], gaia_data['pmdec'][idx_ob2], lw=0, s=4, alpha=1., c='C0')
        # plt.scatter(ves_263['ra'], ves_263['pmdec'], lw=0, s=4, alpha=1., c='C3')
        # plt.xlabel('RA')
        # plt.ylabel('pmDEC')
        # plt.ylim(-7, -2)
        # plt.tight_layout()
        # plt.savefig('aa_pmdec_ra_points_{:.2f}_{:.2f}.png'.format(mean_parallax, d_parallax), dpi=350)
        # plt.close()

        # plt.scatter(gaia_data['dec'][idx_ob2], gaia_data['pmdec'][idx_ob2], lw=0, s=4, alpha=1., c='C0')
        # plt.scatter(ves_263['dec'], ves_263['pmdec'], lw=0, s=4, alpha=1., c='C3')
        ax[1].errorbar(gaia_data['dec'][idx_ob2], gaia_data['pmdec'][idx_ob2], yerr=gaia_data['pmdec_error'][idx_ob2],
                       elinewidth=0.8, marker='.', ms=7, mew=0, alpha=0.75, c='green', ls='None', label='Cyg OB2',
                       capsize=0, capthick=5)
        ax[1].errorbar(ves_263['dec'], ves_263['pmdec'], yerr=ves_263['pmdec_error'],
                       elinewidth=0.8, marker='X', ms=8, mew=0, alpha=0.75, c='red', ls='None', label='VES 263',
                       capsize=0, capthick=5)
        lin_fit = np.polyfit(gaia_data['dec'][idx_ob2], gaia_data['pmdec'][idx_ob2],
                             deg=1, w=1./gaia_data['pmdec_error'][idx_ob2])
        idx_fit_use = np.logical_and(idx_ob2,
                                     np.abs((np.polyval(lin_fit, gaia_data['dec']) - gaia_data['pmdec'])) <= 1.)
        lin_fit = np.polyfit(gaia_data['dec'][idx_fit_use], gaia_data['pmdec'][idx_fit_use],
                             deg=1, w=1. / gaia_data['pmdec_error'][idx_fit_use])
        x_lin = np.linspace(40., 42., 25)
        ax[1].plot(x_lin, np.polyval(lin_fit, x_lin), lw=1, ls='--', alpha=0.75, c='black')
        ax[1].set(xlabel=u'Dec (deg)', ylabel=u'pmDec (mas yr$^{-1}$)', ylim=(-5.7, -2.9), xlim=(40.55, 41.78),
                  xticks=[40.7, 41., 41.3, 41.6], xticklabels=['40.7', '41', '41.3', '41.6'],
                  yticks=[-3, -4, -5], yticklabels=['-3', '-4', '-5'])
        ax[1].grid(ls='--', c='black', alpha=0.2)
        ax[1].tick_params(axis="y", direction="in")
        ax[1].tick_params(axis="x", direction="in")

        plt.tight_layout()
        plt.savefig('linfit_pmra_pmdec_points_{:.2f}_{:.2f}_ver2.pdf'.format(mean_parallax, d_parallax), dpi=300)
        plt.close()

    # plt.scatter(gaia_data['dec'][idx_ob2], gaia_data['pmra'][idx_ob2], lw=0, s=4, alpha=1., c='C0')
    # plt.scatter(ves_263['dec'], ves_263['pmra'], lw=0, s=4, alpha=1., c='C3')
    # plt.xlabel('DEC')
    # plt.ylabel('pmRA')
    # plt.ylim(-6, 0)
    # plt.tight_layout()
    # plt.savefig('aa_pmra_dec_points_{:.2f}_{:.2f}.png'.format(mean_parallax, d_parallax), dpi=350)
    # plt.close()

    # raise SystemExit
    #
    # plt.hist(gaia_data['parallax'], range=(0.45, 0.75), bins=30)
    # plt.hist(gaia_data['parallax'][idx_ob2], range=(0.45, 0.75), bins=30)
    # plt.axvline(ves_263['parallax'], c='blue')
    # plt.show()
    # plt.close()
    #
    # plt.hist(gaia_data['rv'], range=(-60, 40), bins=40)
    # plt.hist(gaia_data['rv'][idx_ob2], range=(-60, 40), bins=40)
    # plt.show()
    # plt.close()
