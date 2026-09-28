// gmatrix - build the sensitivity (Jacobian) matrix G = d(t_pred)/d(model_param)
//           for a model taken from an mcmc_eq output chain, by finite differences
//           through the SAME forward code the inversion uses (forward.c).
//
// The user can then form a model covariance matrix, e.g.
//     Cm = sigma^2 (G^T G)^-1                (unweighted)
//     Cm = (G^T Cd^-1 G)^-1                  (data-weighted, Cd from pick noise)
//
// Version 1: PER-EVENT location G with 4 parameters per event:
//     [ eqx, eqy, eqz, origin_time ]
// Rows are that event's picks (P then S). Derivatives w.r.t. x,y,z are central
// finite differences evaluated through cal_fit_newx/traveltimet, so they are
// exactly consistent with the travel-time tables the model uses. The origin-time
// column is analytic: d(t_pred)/d(ot) = 1 for every pick of the event.
//
// The code is deliberately structured so that a global assembly (all events +
// shared velocity and station-correction columns, "1c") is a straightforward
// extension: see build_event_location_block() and the notes at the bottom.
//
//   SPDX-License-Identifier: GPL-3.0-only
//   Part of the mcmc_eq package.

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "mc.h"

int TRIA, DR, aflag;      // referenced by forward.c / mod_grd.c (globals)
FILE *fpout;              // forward.c references fpout only in print paths we don't call

// ---- prototypes provided by mod_grd.o ----
void read_single_line(FILE *fp, char x[]);
void copy_model(struct Model *dest, struct Model *src);
int  find_in_cell(struct Model *modx, float x);
int  find_neighbor_cell(struct Model *modx, int n);

// ---- prototypes provided by forward.c (included below) ----
// (declared here so we can call them before the include expands)
float cal_fit_newx(struct Model *m, struct DATA *d, float ***tttp, float ***ttts,
                   struct GRDHEAD gh, int calct,
                   float *mfp0, float *mfs0, float *mfp1, float *mfs1,
                   float *mfp2, float *mfs2, float *mfp3, float *mfs3,
                   int flag, int eikonal, int out);
float traveltimet(float **ttt, int nx, int ny, int nz, float h, float dist, float z, float z0);
float dst(float x1, float x2, float y1, float y2);
void  setup_table_new(struct Model *m, float ***ttt, struct GRDHEAD gh, int ps);

// ---- local helpers (also defined in mcmc_eq.c; duplicated here so gmatrix is standalone) ----

// allocate a 3D array of floats [nx][ny][nz]
static float ***make_3d_array(int nx, int ny, int nz)
{
    float ***arr;
    int i, j;
    arr = (float ***) malloc(nx * sizeof(float **));
    for (i = 0; i < nx; i++) {
        arr[i] = (float **) malloc(ny * sizeof(float *));
        for (j = 0; j < ny; j++)
            arr[i][j] = (float *) malloc(nz * sizeof(float));
    }
    return arr;
}

// read the pick file into struct DATA (verbatim behaviour of read_mcmcdata in mcmc_eq.c)
static int read_mcmcdata(FILE *f, struct DATA *d)
{
    int i, ne, end, k;
    char buffer[4096];
    char dum[64], *c;
    int st_id;
    float x, y, z;
    double t;
    char ph[2];
    int np, ns;
    int cl;

    np = 0; ns = 0;
    i = 0; end = 1; k = 0;
    int j = 0;
    while (!feof(f) && end) {
        c = fgets(buffer, 4096, f);
        if (!c) break;
        if (strchr(buffer, '#') == NULL) {
            if (i == MAX_NOQ) { fprintf(stderr, "too many quakes! MAX_NOQ\n"); exit(0); }
            memset(ph, '\0', 2);
            sscanf(buffer, "%s %d %s %f %f %f %lf %d\n",
                   dum, &st_id, ph, &x, &y, &z, &t, &cl);
            if (strchr(ph, 'P') != NULL) {
                if (np == MAX_OBS)   { fprintf(stderr, "too many picks! MAX_OBS\n"); exit(0); }
                if (st_id >= MAX_STAT) { fprintf(stderr, "too many stations! MAX_STAT\n"); exit(0); }
                d[i].p_picks[np].st_id = st_id;
                d[i].p_picks[np].x = x; d[i].p_picks[np].y = y; d[i].p_picks[np].z = z;
                d[i].p_picks[np].t = (float)t; d[i].p_picks[np].cl = cl;
                if (cl == 0) d[i].nobs_p0++; if (cl == 1) d[i].nobs_p1++;
                if (cl == 2) d[i].nobs_p2++; if (cl == 3) d[i].nobs_p3++;
                np++;
            } else {
                if (ns == MAX_OBS)   { fprintf(stderr, "too many picks! MAX_OBS\n"); exit(0); }
                if (st_id >= MAX_STAT) { fprintf(stderr, "too many stations! MAX_STAT\n"); exit(0); }
                d[i].s_picks[ns].st_id = st_id;
                d[i].s_picks[ns].x = x; d[i].s_picks[ns].y = y; d[i].s_picks[ns].z = z;
                d[i].s_picks[ns].t = (float)t; d[i].s_picks[ns].cl = cl;
                if (cl == 0) d[i].nobs_s0++; if (cl == 1) d[i].nobs_s1++;
                if (cl == 2) d[i].nobs_s2++; if (cl == 3) d[i].nobs_s3++;
                ns++;
            }
            k++;
        } else {
            if (j != 0) i++;
            d[i].nobs_p0 = d[i].nobs_p1 = d[i].nobs_p2 = d[i].nobs_p3 = 0;
            d[i].nobs_s0 = d[i].nobs_s1 = d[i].nobs_s2 = d[i].nobs_s3 = 0;
            d[i].xfix = -9999.0; d[i].yfix = -9999.0; d[i].zfix = -9999.0;
            sscanf(buffer, "%s %d %d %d %lf %lf %lf %lf", dum, &d[i].eq_id,
                   &d[i].nobs_p, &d[i].nobs_s, &d[i].reftime, &d[i].xfix, &d[i].yfix, &d[i].zfix);
            k = 0; np = 0; ns = 0;
        }
        j++;
    }
    ne = i + 1;
    return ne;
}

// bring in the exact forward model (cal_fit_newx, traveltimet, setup_table_new, dst)
#include "forward.c"

// ---------------------------------------------------------------------------
// Model reconstruction from a chain line ("bat BF ..." or a chosen "mod ...").
// Line layout (see analyse_eq.c / print_model_raw):
//   <tag> <TY> <number> <dimension> <rms> <p0 p1 p2 p3 s0 s1 s2 s3> then
//   <dimension> triples of (z vp vpvs)
// The EQ and RES lines carry hypocentres and station corrections for that model.
// ---------------------------------------------------------------------------

// parse the model (velocity + noise + dimension) from a "bat"/"mod" token line
static int parse_model_line(char *buf, struct Model *m)
{
    char *tok;
    tok = strtok(buf, " ");            // tag: bat/mod/sta
    if (!tok) return 0;
    strtok(NULL, " ");                 // TY
    m->number    = atol(strtok(NULL, " "));
    m->dimension = atol(strtok(NULL, " "));
    strtok(NULL, " ");                 // rms
    m->p_noise0 = atof(strtok(NULL, " ")); m->p_noise1 = atof(strtok(NULL, " "));
    m->p_noise2 = atof(strtok(NULL, " ")); m->p_noise3 = atof(strtok(NULL, " "));
    m->s_noise0 = atof(strtok(NULL, " ")); m->s_noise1 = atof(strtok(NULL, " "));
    m->s_noise2 = atof(strtok(NULL, " ")); m->s_noise3 = atof(strtok(NULL, " "));
    for (long l = 0; l < m->dimension; l++) {
        m->z[l]    = atof(strtok(NULL, " "));
        m->vp[l]   = atof(strtok(NULL, " "));
        m->vpvs[l] = atof(strtok(NULL, " "));
    }
    return 1;
}

// Read ONE model set from the chain identified by tag ("bat" for best-fit, or the
// last "mod" line if tag=="mod"). Fills m (velocity+noise), eq hypocentres, and
// station corrections (pres/sres). We use the model whose "number" matches the
// chosen model line so EQ/RES belong to the same sample.
static int load_model_set(const char *chain, const char *want_tag,
                          struct Model *m, int *noq_out, int *nos_out)
{
    FILE *f = fopen(chain, "r");
    if (!f) { fprintf(stderr, "cannot open chain %s\n", chain); return 0; }

    char buf[200000];
    char keep_model[200000];
    long target_number = -1;
    int have_model = 0;
    int noq = 0, nos = 0;

    // init station corrections to "unset"
    for (int i = 0; i < MAX_STAT; i++) { m->pres[i] = -99999; m->sres[i] = -99999; }

    // Pass 1: find the desired model line (its dimension/velocity/noise) and its "number".
    while (fgets(buf, sizeof(buf), f)) {
        if (strncmp(buf, want_tag, strlen(want_tag)) == 0) {
            strcpy(keep_model, buf);
            // peek number without destroying: quick parse
            char tmp[200000]; strcpy(tmp, buf);
            strtok(tmp, " "); strtok(NULL, " ");
            target_number = atol(strtok(NULL, " "));
            have_model = 1;
            if (strcmp(want_tag, "bat") == 0) break;  // best-fit is unique
            // for "mod" we keep the LAST occurrence
        }
    }
    if (!have_model) { fprintf(stderr, "tag '%s' not found in chain\n", want_tag); fclose(f); return 0; }
    parse_model_line(keep_model, m);

    // Pass 2: collect EQ and RES lines belonging to target_number.
    rewind(f);
    while (fgets(buf, sizeof(buf), f)) {
        if (strncmp(buf, "EQ", 2) == 0) {
            char tmp[4096]; strncpy(tmp, buf, sizeof(tmp)-1); tmp[sizeof(tmp)-1]=0;
            strtok(tmp, " "); strtok(NULL, " ");
            long num = atol(strtok(NULL, " "));
            int  eqi = atoi(strtok(NULL, " "));
            strtok(NULL, " ");                          // rms
            float ex = atof(strtok(NULL, " "));
            float ey = atof(strtok(NULL, " "));
            float ez = atof(strtok(NULL, " "));
            if (num == target_number) {
                m->eq[eqi].x = ex; m->eq[eqi].y = ey; m->eq[eqi].z = ez;
                if (eqi + 1 > noq) noq = eqi + 1;
            }
        } else if (strncmp(buf, "RES", 3) == 0) {
            char tmp[4096]; strncpy(tmp, buf, sizeof(tmp)-1); tmp[sizeof(tmp)-1]=0;
            strtok(tmp, " "); strtok(NULL, " ");
            long num = atol(strtok(NULL, " "));
            int  sti = atoi(strtok(NULL, " "));
            strtok(NULL, " ");                          // rms
            float pr = atof(strtok(NULL, " "));
            float sr = atof(strtok(NULL, " "));
            if (num == target_number) {
                m->pres[sti] = pr; m->sres[sti] = sr;
                if (sti + 1 > nos) nos = sti + 1;
            }
        }
    }
    fclose(f);
    *noq_out = noq; *nos_out = nos;
    return 1;
}

// ---------------------------------------------------------------------------
// Predicted travel time for a single pick at an arbitrary hypocentre (ex,ey,ez),
// using the exact forward expressions from cal_fit_newx. This lets us perturb
// one hypocentre coordinate and recompute t_pred with identical arithmetic,
// without re-running the whole ensemble loop.
// tables tttp/ttts must be current for the model m.
// NOTE: does NOT include the origin-time term (that is handled analytically as a
// separate G column). It returns TT(dist,ez) + station_correction, matching the
// t_pred stored by cal_fit_newx.
// ---------------------------------------------------------------------------
static float pred_time_pick(struct Model *m, struct OBS *pk, int is_s,
                            float ex, float ey, float ez,
                            float ***tttp, float ***ttts, struct GRDHEAD gh,
                            int eikonal)
{
    float dist = dst(pk->x, ex, pk->y, ey);
    float tt;
    if (eikonal == 0) {
        int sc = find_in_cell(m, 0.0);
        float v = is_s ? (m->vp[sc] / m->vpvs[sc]) : m->vp[sc];
        tt = sqrtf(dist * dist + ez * ez) / v;
    } else {
        float ***T = is_s ? ttts : tttp;
        tt = traveltimet(T[pk->layer],   gh.nx, gh.ny, gh.nz, gh.h, dist, ez, gh.z0) * pk->w1
           + traveltimet(T[pk->layer+1], gh.nx, gh.ny, gh.nz, gh.h, dist, ez, gh.z0) * pk->w2;
    }
    float cor = is_s ? m->sres[pk->st_id] : m->pres[pk->st_id];
    return tt + cor;
}

// ---------------------------------------------------------------------------
// Per-event location G block.
// For event i: rows = its P picks then S picks; columns = [d/dx, d/dy, d/dz, d/dot].
//   x,y,z columns: central finite difference through pred_time_pick.
//   ot  column   : analytic, = +1 for every pick (t_pred increases 1:1 with ot).
// The block is written to fp with a header identifying the event and each row's
// station id / phase / observed noise (so weights can be applied downstream).
// Returns number of rows written.
// ---------------------------------------------------------------------------
static int build_event_location_block(FILE *fp, struct Model *m, struct DATA *d,
                                       int evt, float ***tttp, float ***ttts,
                                       struct GRDHEAD gh, int eikonal,
                                       float dx, float dz)
{
    float ex = m->eq[evt].x, ey = m->eq[evt].y, ez = m->eq[evt].z;
    int npk = d[evt].nobs_p + d[evt].nobs_s;
    fprintf(fp, "# EVENT %d  nrows %d  ncols 4  params[x y z ot]  hypo %.4f %.4f %.4f\n",
            evt, npk, ex, ey, ez);

    // per-pick noise: choose the class-appropriate noise (sigma) so the user can
    // build Cd^-1 = diag(1/sigma^2). We report it in the row header.
    float pn[4] = { m->p_noise0, m->p_noise1, m->p_noise2, m->p_noise3 };
    float sn[4] = { m->s_noise0, m->s_noise1, m->s_noise2, m->s_noise3 };

    for (int phase = 0; phase < 2; phase++) {          // 0 = P, 1 = S
        int nobs = phase ? d[evt].nobs_s : d[evt].nobs_p;
        struct OBS *picks = phase ? d[evt].s_picks : d[evt].p_picks;
        for (int j = 0; j < nobs; j++) {
            struct OBS *pk = &picks[j];
            // central differences
            float gx = (pred_time_pick(m, pk, phase, ex + dx, ey, ez, tttp, ttts, gh, eikonal)
                      - pred_time_pick(m, pk, phase, ex - dx, ey, ez, tttp, ttts, gh, eikonal)) / (2.0f * dx);
            float gy = (pred_time_pick(m, pk, phase, ex, ey + dx, ez, tttp, ttts, gh, eikonal)
                      - pred_time_pick(m, pk, phase, ex, ey - dx, ez, tttp, ttts, gh, eikonal)) / (2.0f * dx);
            float gz = (pred_time_pick(m, pk, phase, ex, ey, ez + dz, tttp, ttts, gh, eikonal)
                      - pred_time_pick(m, pk, phase, ex, ey, ez - dz, tttp, ttts, gh, eikonal)) / (2.0f * dz);
            float got = 1.0f;                          // analytic origin-time derivative
            float sigma = phase ? sn[pk->cl] : pn[pk->cl];
            fprintf(fp, "%d %c st %d cl %d sigma %.5f  G %.8e %.8e %.8e %.8e\n",
                    evt, phase ? 'S' : 'P', pk->st_id, pk->cl, sigma, gx, gy, gz, got);
        }
    }
    return npk;
}

// predicted time for a pick at the model's stored hypocentre, using given tables.
// (thin wrapper used by velocity-perturbation columns: hypocentre is fixed, only
//  the tables/velocities change.)
static float pred_time_pick_base(struct Model *m, struct OBS *pk, int is_s,
                                 int evt, float ***tttp, float ***ttts,
                                 struct GRDHEAD gh, int eikonal)
{
    return pred_time_pick(m, pk, is_s, m->eq[evt].x, m->eq[evt].y, m->eq[evt].z,
                          tttp, ttts, gh, eikonal);
}

// ---------------------------------------------------------------------------
// GLOBAL 1c assembly.
// Parameter vector (columns), in this fixed order:
//   loc:  4*ne         -> [x_e, y_e, z_e, ot_e] for e = 0..ne-1
//   vp:   dim          -> vp_k     for k = 0..dim-1
//   vpvs: dim          -> vpvs_k   for k = 0..dim-1
//   pres: nos          -> pres_s   for s = 0..nos-1   (P station corrections)
//   sres: nos          -> sres_s   for s = 0..nos-1   (S station corrections)
// Rows: every P then S pick of every event (data rows), followed by any
//   constraint pseudo-rows (station-correction zero-mean, if scor_flag==0).
//
// Output format (single file):
//   # manifest header lines (column ranges, dims, sigmas, constraint info)
//   D <row> <phase> <evt> <st> <cl> <sigma> <resid>      (one per data row)
//   G <row> <col> <value>                                (sparse triplets)
//   C <row> <kind> <weight>                              (constraint row meta)
// Downstream you assemble the sparse G, weight matrix W=diag(1/sigma^2) (data)
// or the constraint weight (pseudo-rows), and form G^T W G etc.
//
// Column derivatives:
//   x,y,z : central FD via pred_time_pick (hypocentre perturbed), block-diagonal.
//   ot    : analytic +1 for that event's picks.
//   pres_s: analytic +1 for a P pick recorded at station s (0 else).
//   sres_s: analytic +1 for an S pick recorded at station s (0 else).
//   vp_k  : central FD; perturb m.vp[k], rebuild BOTH tables, recompute all t_pred.
//   vpvs_k: central FD; perturb m.vpvs[k], rebuild S table only.
// ---------------------------------------------------------------------------
static int build_global_G(FILE *fp, struct Model *m, struct DATA *d, int ne, int nos,
                          float ***tttp, float ***ttts, struct GRDHEAD gh,
                          int eikonal, float dx, float dz, float dvp, float dvpvs,
                          int scor_flag, int reference_station, float cweight,
                          int fix_vp, int fix_vpvs, int fix_loc, int fix_res)
{
    int dim = (int)m->dimension;
    // Per-event location parameters: normally [x y z ot] (4). If quake locations
    // were held fixed, x/y/z are not free -- keep only the origin-time column
    // (1/event), which mcmc_eq still sets analytically per event.
    int loc_per = fix_loc ? 1 : 4;
    int nvel    = fix_vp   ? 0 : dim;    // vp columns   (0 if vp not inverted)
    int nvpvs   = fix_vpvs ? 0 : dim;    // vpvs columns (0 if vp/vs not inverted)
    int npres   = fix_res ? 0 : nos;     // station P-correction columns
    int nsres   = fix_res ? 0 : nos;     // station S-correction columns

    // column offsets (blocks with zero width simply collapse)
    int c_loc  = 0;                      // loc_per*ne location columns
    int c_vp   = c_loc  + loc_per * ne;  // nvel  vp columns
    int c_vpvs = c_vp   + nvel;          // nvpvs vpvs columns
    int c_pres = c_vpvs + nvpvs;         // npres pres columns
    int c_sres = c_pres + npres;         // nsres sres columns
    int ncols  = c_sres + nsres;

    float pn[4] = { m->p_noise0, m->p_noise1, m->p_noise2, m->p_noise3 };
    float sn[4] = { m->s_noise0, m->s_noise1, m->s_noise2, m->s_noise3 };

    // ---- assign a global row index to every pick, and remember its identity ----
    // We also store its predicted time (base) so we can difference for velocity cols.
    int nrows_data = 0;
    for (int e = 0; e < ne; e++) nrows_data += d[e].nobs_p + d[e].nobs_s;

    // per-row bookkeeping
    int   *r_evt   = malloc(nrows_data * sizeof(int));
    int   *r_phase = malloc(nrows_data * sizeof(int));   // 0=P,1=S
    int   *r_pkidx = malloc(nrows_data * sizeof(int));   // index within that phase's pick array
    float *tbase   = malloc(nrows_data * sizeof(float)); // base predicted time
    if (!r_evt || !r_phase || !r_pkidx || !tbase) { fprintf(stderr,"malloc rows failed\n"); return -1; }

    int row = 0;
    for (int e = 0; e < ne; e++) {
        for (int ph = 0; ph < 2; ph++) {
            int nobs = ph ? d[e].nobs_s : d[e].nobs_p;
            struct OBS *pk = ph ? d[e].s_picks : d[e].p_picks;
            for (int j = 0; j < nobs; j++) {
                r_evt[row] = e; r_phase[row] = ph; r_pkidx[row] = j;
                tbase[row] = pred_time_pick_base(m, &pk[j], ph, e, tttp, ttts, gh, eikonal);
                row++;
            }
        }
    }

    // ---- manifest header ----
    fprintf(fp, "# gmatrix v2 GLOBAL 1c  ncols %d  ndata_rows %d\n", ncols, nrows_data);
    // colmap: only list blocks that actually have columns. loc uses loc_per/event
    // (4 = x y z ot, or 1 = ot only when quake locations are fixed).
    fprintf(fp, "# colmap loc %d..%d (%d/event: %s)",
            c_loc, c_loc + loc_per*ne - 1, loc_per, fix_loc ? "ot" : "x y z ot");
    if (nvel)  fprintf(fp, "  vp %d..%d", c_vp, c_vp+nvel-1);
    if (nvpvs) fprintf(fp, "  vpvs %d..%d", c_vpvs, c_vpvs+nvpvs-1);
    if (npres) fprintf(fp, "  pres %d..%d  sres %d..%d", c_pres, c_pres+npres-1, c_sres, c_sres+nsres-1);
    fprintf(fp, "\n");
    fprintf(fp, "# dims ne %d  dim %d  nos %d  eikonal %d  TRIA %d\n", ne, dim, nos, eikonal, TRIA);
    fprintf(fp, "# fixed  vp %d  vpvs %d  loc %d  res %d   (1 = held fixed / not a free parameter)\n",
            fix_vp, fix_vpvs, fix_loc, fix_res);
    fprintf(fp, "# fd dx %.4f dz %.4f dvp %.4f dvpvs %.4f  scor_flag %d ref_station %d cweight %.3e\n",
            dx, dz, dvp, dvpvs, scor_flag, reference_station, cweight);
    fprintf(fp, "# rows: D=data row meta, G=sparse triplet (row col value), C=constraint row meta\n");

    // ---- data-row meta (residual + sigma) ----
    for (int r = 0; r < nrows_data; r++) {
        int e = r_evt[r], ph = r_phase[r], j = r_pkidx[r];
        struct OBS *pk = ph ? &d[e].s_picks[j] : &d[e].p_picks[j];
        float sigma = ph ? sn[pk->cl] : pn[pk->cl];
        // residual = observed - predicted - origin_time_shift (event mean already in m->origin)
        float resid = pk->t - tbase[r] + m->origin[e];  // t_pred includes station corr; origin from cal_fit
        fprintf(fp, "D %d %c %d %d %d %.6f %.6f\n", r, ph ? 'S' : 'P', e, pk->st_id, pk->cl, sigma, resid);
    }

    // ---- LOCATION columns (central FD) + analytic ot + analytic station cols ----
    for (int r = 0; r < nrows_data; r++) {
        int e = r_evt[r], ph = r_phase[r], j = r_pkidx[r];
        struct OBS *pk = ph ? &d[e].s_picks[j] : &d[e].p_picks[j];
        int base = c_loc + loc_per*e;
        if (!fix_loc) {
            // free hypocentre: x,y,z via central FD, then analytic origin-time
            float ex = m->eq[e].x, ey = m->eq[e].y, ez = m->eq[e].z;
            float gx = (pred_time_pick(m, pk, ph, ex+dx, ey, ez, tttp, ttts, gh, eikonal)
                      - pred_time_pick(m, pk, ph, ex-dx, ey, ez, tttp, ttts, gh, eikonal)) / (2.0f*dx);
            float gy = (pred_time_pick(m, pk, ph, ex, ey+dx, ez, tttp, ttts, gh, eikonal)
                      - pred_time_pick(m, pk, ph, ex, ey-dx, ez, tttp, ttts, gh, eikonal)) / (2.0f*dx);
            float gz = (pred_time_pick(m, pk, ph, ex, ey, ez+dz, tttp, ttts, gh, eikonal)
                      - pred_time_pick(m, pk, ph, ex, ey, ez-dz, tttp, ttts, gh, eikonal)) / (2.0f*dz);
            if (gx != 0.0f) fprintf(fp, "G %d %d %.8e\n", r, base+0, gx);
            if (gy != 0.0f) fprintf(fp, "G %d %d %.8e\n", r, base+1, gy);
            if (gz != 0.0f) fprintf(fp, "G %d %d %.8e\n", r, base+2, gz);
            fprintf(fp, "G %d %d %.8e\n", r, base+3, 1.0);  // d/d origin_time
        } else {
            // quake locations fixed: only the analytic origin-time column remains
            fprintf(fp, "G %d %d %.8e\n", r, base+0, 1.0);  // d/d origin_time
        }

        // station-correction column (analytic +1), unless corrections are fixed or
        // that phase/station is frozen by scor_flag
        if (!fix_res) {
            int st = pk->st_id;
            if (ph == 0) {  // P pick -> pres[st]
                int drop = (scor_flag == 1 || scor_flag == 2) && (st == reference_station);
                int omit = (scor_flag == -2);   // S-only inversion: P corrections not free
                if (!drop && !omit && st < nos) fprintf(fp, "G %d %d %.8e\n", r, c_pres + st, 1.0);
            } else {        // S pick -> sres[st]
                int drop = (scor_flag == 2) && (st == reference_station);
                int omit = (scor_flag == -1);   // P-only inversion: S corrections not free
                if (!drop && !omit && st < nos) fprintf(fp, "G %d %d %.8e\n", r, c_sres + st, 1.0);
            }
        }
    }

    // ---- VELOCITY columns: perturb model param, rebuild tables, recompute t_pred ----
    // vp and vp/vs are gated independently: vp columns are skipped when vp was not
    // inverted (no P/M/B/D proposals), vp/vs columns when vp/vs was not inverted
    // (no 'V' proposal). Central difference throughout; tttp/ttts hold the BASE
    // tables, perturbed tables go to pP/pS scratch so the base stays valid.
    if ((!fix_vp || !fix_vpvs) && eikonal == 1 && dim > 0) {
        int nxmod = (int)sqrtf((float)gh.nx*gh.nx + (float)gh.ny*gh.ny);
        float ***pP = make_3d_array(gh.nz, gh.nz, nxmod);
        float ***pS = make_3d_array(gh.nz, gh.nz, nxmod);
        float *tplus  = malloc(nrows_data * sizeof(float));
        float *tminus = malloc(nrows_data * sizeof(float));
        if (!pP || !pS || !tplus || !tminus) { fprintf(stderr,"malloc vel failed\n"); return -1; }

        // vp_k : vs = vp/vpvs, so a change in vp affects BOTH P and S -> rebuild both.
        for (int k = 0; !fix_vp && k < dim; k++) {
            float v0 = m->vp[k];

            m->vp[k] = v0 + dvp;
            setup_table_new(m, pP, gh, 1);
            setup_table_new(m, pS, gh, 2);
            for (int r = 0; r < nrows_data; r++) {
                int e = r_evt[r], ph = r_phase[r], j = r_pkidx[r];
                struct OBS *pk = ph ? &d[e].s_picks[j] : &d[e].p_picks[j];
                tplus[r] = pred_time_pick_base(m, pk, ph, e, pP, pS, gh, eikonal);
            }
            m->vp[k] = v0 - dvp;
            setup_table_new(m, pP, gh, 1);
            setup_table_new(m, pS, gh, 2);
            for (int r = 0; r < nrows_data; r++) {
                int e = r_evt[r], ph = r_phase[r], j = r_pkidx[r];
                struct OBS *pk = ph ? &d[e].s_picks[j] : &d[e].p_picks[j];
                tminus[r] = pred_time_pick_base(m, pk, ph, e, pP, pS, gh, eikonal);
            }
            m->vp[k] = v0;   // restore

            int col = c_vp + k;
            for (int r = 0; r < nrows_data; r++) {
                float g = (tplus[r] - tminus[r]) / (2.0f * dvp);
                if (g != 0.0f) fprintf(fp, "G %d %d %.8e\n", r, col, g);
            }
            fprintf(stderr, "\rvelocity col vp_%d/%d", k+1, dim);
        }

        // vpvs_k : affects S velocity only (vs = vp/vpvs) -> rebuild S table only.
        for (int k = 0; !fix_vpvs && k < dim; k++) {
            float v0 = m->vpvs[k];

            m->vpvs[k] = v0 + dvpvs;
            setup_table_new(m, pS, gh, 2);
            for (int r = 0; r < nrows_data; r++) {
                int e = r_evt[r], ph = r_phase[r], j = r_pkidx[r];
                if (ph == 0) { tplus[r] = 0.0f; continue; }  // P picks insensitive to vpvs
                struct OBS *pk = &d[e].s_picks[j];
                tplus[r] = pred_time_pick_base(m, pk, 1, e, tttp, pS, gh, eikonal);
            }
            m->vpvs[k] = v0 - dvpvs;
            setup_table_new(m, pS, gh, 2);
            for (int r = 0; r < nrows_data; r++) {
                int e = r_evt[r], ph = r_phase[r], j = r_pkidx[r];
                if (ph == 0) { tminus[r] = 0.0f; continue; }
                struct OBS *pk = &d[e].s_picks[j];
                tminus[r] = pred_time_pick_base(m, pk, 1, e, tttp, pS, gh, eikonal);
            }
            m->vpvs[k] = v0;   // restore

            int col = c_vpvs + k;
            for (int r = 0; r < nrows_data; r++) {
                if (r_phase[r] == 0) continue;  // exact zero for P
                float g = (tplus[r] - tminus[r]) / (2.0f * dvpvs);
                if (g != 0.0f) fprintf(fp, "G %d %d %.8e\n", r, col, g);
            }
            fprintf(stderr, "\rvelocity col vpvs_%d/%d", k+1, dim);
        }
        fprintf(stderr, "\n");

        // base tables (tttp/ttts) were never modified; nothing to restore.
        for (int a=0;a<gh.nz;a++){for(int b=0;b<gh.nz;b++){free(pP[a][b]);free(pS[a][b]);}free(pP[a]);free(pS[a]);}
        free(pP); free(pS); free(tplus); free(tminus);
    }

    // ---- constraint pseudo-rows for zero-mean station corrections (scor_flag==0) ----
    // Only meaningful when station-correction columns exist (not fixed).
    int extra = 0;
    if (scor_flag == 0 && !fix_res) {
        int rP = nrows_data + extra;
        fprintf(fp, "C %d meanP %.6e\n", rP, cweight);
        for (int s = 0; s < nos; s++) fprintf(fp, "G %d %d %.8e\n", rP, c_pres + s, cweight);
        extra++;
        int rS = nrows_data + extra;
        fprintf(fp, "C %d meanS %.6e\n", rS, cweight);
        for (int s = 0; s < nos; s++) fprintf(fp, "G %d %d %.8e\n", rS, c_sres + s, cweight);
        extra++;
    }

    free(r_evt); free(r_phase); free(r_pkidx); free(tbase);
    return nrows_data + extra;
}

// ---------------------------------------------------------------------------
// Read only the fields we need from the config: grid header + eikonal flag + TRIA
// + the two lines that determine which parameters were actually FREE in the run,
// so we can omit fixed classes' columns from G.
//
// Two independent mechanisms fix parameters:
//  (1) line 33 "model modification tests" (normal MCMC): a letter string listing
//      which proposal moves are enabled. Relevant letters:
//        'V'          -> vp/vs perturbation (S table)            => vpvs is FREE
//        'P','M','B','D' -> vp / node moves (P & S tables)       => vp is FREE
//        'Q' -> quake moves (locations FREE); 'R' -> corrections FREE; 'N' -> noise FREE
//  (2) line 34 "aflag switch" with aflag==3 (JHD/input-from-file): letters in the
//      switch name classes READ FROM FILE and thus held fixed: 'V' vel(vp&vpvs),
//      'Q' locations, 'R' corrections, 'N' noise.
//
// A class is fixed if it is disabled in (1) OR frozen by (2). fix_* outputs
// (1 = fixed / not a free parameter):
//   fix_vp, fix_vpvs : velocity components (separately!)
//   fix_loc          : quake locations (x,y,z)
//   fix_res          : station corrections
// ---------------------------------------------------------------------------
static int read_config(const char *path, struct GRDHEAD *gh, int *eikonal, int *max_dim,
                       int *scor_flag, int *reference_station,
                       int *fix_vp, int *fix_vpvs, int *fix_loc, int *fix_res)
{
    FILE *cf = fopen(path, "r");
    if (!cf) { fprintf(stderr, "cannot open config %s\n", path); return 0; }
    char line[1024];
    int idum; int tr;
    float rp, rs;
    read_single_line(cf, line); sscanf(line, "%f", &gh->h);
    read_single_line(cf, line); sscanf(line, "%d", &gh->nx);
    read_single_line(cf, line); sscanf(line, "%d", &gh->ny);
    read_single_line(cf, line); sscanf(line, "%d", &gh->nz);
    read_single_line(cf, line); sscanf(line, "%f", &gh->x0);
    read_single_line(cf, line); sscanf(line, "%f", &gh->y0);
    read_single_line(cf, line); sscanf(line, "%f", &gh->z0);
    read_single_line(cf, line); sscanf(line, "%d", max_dim);
    // 19 scalar lines: vpmin,vpmax,vpvsmin,vpvsmax,noise_min,noise_max,res_min,
    // res_max,sdevx,sdevy,sdevz,sdevvp,sdevvpvs,sdevn,(sdevxs epi),sdevys,sdevzs,
    // sdevresidual,inv_control
    for (int i = 0; i < 19; i++) read_single_line(cf, line);
    // reference_station line: "reference_station scor_flag ref_statcor_P ref_statcor_S"
    read_single_line(cf, line); sscanf(line, "%d %d %f %f", reference_station, scor_flag, &rp, &rs);
    read_single_line(cf, line); sscanf(line, "%d", &tr); TRIA = tr;   // TRIA (Voronoi=0)
    read_single_line(cf, line);                                       // number of models in chain
    read_single_line(cf, line);                                       // deci (output every nth)
    read_single_line(cf, line); sscanf(line, "%d %d", &idum, eikonal);// true_random eikonal

    // ---- determine which parameters were free ----
    // line 33: model modification tests (normal-MCMC proposal letters)
    // line 34: aflag + JHD switch (input-from-file fixes classes when aflag==3)
    *fix_vp = *fix_vpvs = *fix_loc = *fix_res = 0;
    char tests[128] = "";
    read_single_line(cf, line);                                       // line 33: tests, e.g. "QN QVRPBDMN"
    strncpy(tests, line, sizeof(tests) - 1);
    read_single_line(cf, line);                                       // line 34: mode, e.g. "0 VRN" or "3 QV"
    {
        int aflag = -1; char sw[64] = "";
        sscanf(line, "%d %63s", &aflag, sw);

        // (1) normal-MCMC freedom from the line-33 test string:
        //     vpvs free iff 'V' present; vp free iff any of P/M/B/D present;
        //     locations free iff 'Q'; corrections free iff 'R'.
        int vpvs_free = (strchr(tests, 'V') != NULL);
        int vp_free   = (strchr(tests, 'P') || strchr(tests, 'M') ||
                         strchr(tests, 'B') || strchr(tests, 'D'));
        int loc_free  = (strchr(tests, 'Q') != NULL);
        int res_free  = (strchr(tests, 'R') != NULL);
        if (!vpvs_free) *fix_vpvs = 1;
        if (!vp_free)   *fix_vp   = 1;
        if (!loc_free)  *fix_loc  = 1;
        if (!res_free)  *fix_res  = 1;

        // (2) JHD / input-from-file (aflag==3): listed classes are read from file
        //     and held fixed, overriding (1) toward "fixed".
        if (aflag == 3) {
            if (strchr(sw, 'V')) { *fix_vp = 1; *fix_vpvs = 1; }  // both vp & vpvs
            if (strchr(sw, 'Q'))   *fix_loc = 1;
            if (strchr(sw, 'R'))   *fix_res = 1;
        }
    }
    fclose(cf);
    return 1;
}

int main(int argc, char *argv[])
{
    if (argc < 4) {
        fprintf(stderr,
          "usage: %s config.dat chain_file pick_file [mode] [tag] [outfile] [fdstep] [cweight]\n"
          "  mode    : 'loc' per-event location G (default) | 'global' full 1c matrix\n"
          "  tag     : 'bat' (best-fit model, default) or 'mod' (last sampled model)\n"
          "  outfile : where to write G (default: gmatrix.out)\n"
          "  fdstep  : finite-difference step in km for x,y,z (default 0.01)\n"
          "  cweight : constraint pseudo-row weight for zero-mean station corr (default 1e6)\n"
          "\n"
          "loc:    per-event 4 params [x y z origin_time]; Cm = sigma^2 (G^T G)^-1 per event.\n"
          "global: params [(x y z ot)/event | vp_k | vpvs_k | pres_s | sres_s].\n"
          "        Emits sparse triplets + data/constraint meta for Cm and resolution R.\n", argv[0]);
        return 1;
    }
    const char *cfg   = argv[1];
    const char *chain = argv[2];
    const char *picks = argv[3];
    const char *mode  = (argc > 4) ? argv[4] : "loc";
    const char *tag   = (argc > 5) ? argv[5] : "bat";
    // default output name embeds the chain id (e.g. rjx-009_4621867.out -> gmatrix_009.out)
    // so matrices from different chains don't overwrite each other. An explicit 6th
    // argument overrides this.
    char outf_buf[512];
    if (argc > 6) {
        snprintf(outf_buf, sizeof(outf_buf), "%s", argv[6]);
    } else {
        // extract chain id: basename, strip "rjx-" prefix, take up to first '_' or '.'
        const char *bn = strrchr(chain, '/');
        bn = bn ? bn + 1 : chain;
        char cid[128] = "";
        const char *p = bn;
        if (strncmp(p, "rjx-", 4) == 0) p += 4;
        int k = 0;
        while (*p && *p != '_' && *p != '.' && k < (int)sizeof(cid) - 1) cid[k++] = *p++;
        cid[k] = '\0';
        if (k > 0) snprintf(outf_buf, sizeof(outf_buf), "gmatrix_%s.out", cid);
        else       snprintf(outf_buf, sizeof(outf_buf), "gmatrix.out");
    }
    const char *outf = outf_buf;
    float fd_step     = (argc > 7) ? atof(argv[7]) : 0.01f;  // FD step (km) for x,y,z
    float cweight     = (argc > 8) ? atof(argv[8]) : 1.0e6f; // zero-mean constraint weight
    float dvp_step    = (argc > 9) ? atof(argv[9]) : 0.02f;  // FD step (km/s) for vp
    float dvpvs_step  = (argc > 10) ? atof(argv[10]) : 0.01f;// FD step for vp/vs
    int   global_mode = (strcmp(mode, "global") == 0);

    struct GRDHEAD gh;
    int eikonal = 1, max_dim = 0, scor_flag = 0, reference_station = 0;
    int fix_vp = 0, fix_vpvs = 0, fix_loc = 0, fix_res = 0;
    if (!read_config(cfg, &gh, &eikonal, &max_dim, &scor_flag, &reference_station,
                     &fix_vp, &fix_vpvs, &fix_loc, &fix_res)) return 1;
    fprintf(stderr, "grid: h=%.3f nx=%d ny=%d nz=%d z0=%.3f  eikonal=%d TRIA=%d scor_flag=%d ref=%d\n",
            gh.h, gh.nx, gh.ny, gh.nz, gh.z0, eikonal, TRIA, scor_flag, reference_station);
    if (fix_vp || fix_vpvs || fix_loc || fix_res)
        fprintf(stderr, "fixed parameters (omitted from G): vp=%d vpvs=%d loc=%d res=%d\n",
                fix_vp, fix_vpvs, fix_loc, fix_res);

    // load picks
    struct DATA *s = malloc(MAX_NOQ * sizeof(struct DATA));
    if (!s) { fprintf(stderr, "malloc s failed\n"); return 1; }
    FILE *pf = fopen(picks, "r");
    if (!pf) { fprintf(stderr, "cannot open pick file %s\n", picks); return 1; }
    int ne = read_mcmcdata(pf, s);
    fclose(pf);
    fprintf(stderr, "read %d events from %s\n", ne, picks);

    // precompute receiver-elevation interpolation weights (identical to mcmc_eq main)
    for (int i = 0; i < ne; i++) {
        for (int j = 0; j < s[i].nobs_p; j++) {
            s[i].p_picks[j].layer = (int)((s[i].p_picks[j].z - gh.z0) / gh.h);
            s[i].p_picks[j].w2 = -(s[i].p_picks[j].layer * gh.h + gh.z0 - s[i].p_picks[j].z) / gh.h;
            s[i].p_picks[j].w1 = 1.0 - s[i].p_picks[j].w2;
        }
        for (int j = 0; j < s[i].nobs_s; j++) {
            s[i].s_picks[j].layer = (int)((s[i].s_picks[j].z - gh.z0) / gh.h);
            s[i].s_picks[j].w2 = -(s[i].s_picks[j].layer * gh.h + gh.z0 - s[i].s_picks[j].z) / gh.h;
            s[i].s_picks[j].w1 = 1.0 - s[i].s_picks[j].w2;
        }
    }

    // reconstruct the chosen model (velocity + hypocentres + station corrections)
    struct Model m;
    memset(&m, 0, sizeof(m));
    int noq = 0, nos = 0;
    if (!load_model_set(chain, tag, &m, &noq, &nos)) return 1;
    m.noq = ne;                     // hypocentre count follows the pick file
    fprintf(stderr, "model '%s' #%ld dim=%ld  events=%d stations=%d\n",
            tag, m.number, m.dimension, noq, nos);
    if (noq < ne)
        fprintf(stderr, "WARNING: chain provided %d EQ lines but pick file has %d events\n", noq, ne);

    // allocate + populate travel-time tables for this model (eikonal path)
    int nxmod = (int) sqrtf((float)gh.nx * gh.nx + (float)gh.ny * gh.ny);
    float ***tttp = make_3d_array(gh.nz, gh.nz, nxmod);
    float ***ttts = make_3d_array(gh.nz, gh.nz, nxmod);
    if (eikonal == 1) {
        setup_table_new(&m, tttp, gh, 1);
        setup_table_new(&m, ttts, gh, 2);
    }

    // sanity: run the full forward once so t_pred is populated and RMS is sensible
    float mfp0, mfs0, mfp1, mfs1, mfp2, mfs2, mfp3, mfs3;
    cal_fit_newx(&m, s, tttp, ttts, gh, 0,
                 &mfp0, &mfs0, &mfp1, &mfs1, &mfp2, &mfs2, &mfp3, &mfs3, 0, eikonal, 0);

    // finite-difference steps: small vs. the grid spacing. Central differences.
    float dx = fd_step;      // km, horizontal & used for x,y
    float dz = fd_step;      // km, depth
    float dvp = dvp_step;    // km/s, vp perturbation
    float dvpvs = dvpvs_step;// dimensionless, vp/vs perturbation

    FILE *fp = fopen(outf, "w");
    if (!fp) { fprintf(stderr, "cannot open output %s\n", outf); return 1; }

    if (global_mode) {
        int total = build_global_G(fp, &m, s, ne, nos, tttp, ttts, gh, eikonal,
                                   dx, dz, dvp, dvpvs, scor_flag, reference_station, cweight,
                                   fix_vp, fix_vpvs, fix_loc, fix_res);
        fclose(fp);
        if (total < 0) { fprintf(stderr, "global assembly failed\n"); return 1; }
        fprintf(stderr, "wrote global G: %d total rows (incl. constraints) to %s\n", total, outf);
    } else {
        if (fix_loc)
            fprintf(stderr, "WARNING: quake locations were held fixed (JHD 'Q'); the per-event "
                            "location covariance from loc-mode is not meaningful for this run.\n");
        fprintf(fp, "# gmatrix v1 per-event location G  cols=[dTp/dx dTp/dy dTp/dz dTp/dot]\n");
        fprintf(fp, "# model tag=%s number=%ld dim=%ld  fd dx=%.4f dz=%.4f eikonal=%d\n",
                tag, m.number, m.dimension, dx, dz, eikonal);
        int total = 0;
        for (int i = 0; i < ne; i++)
            total += build_event_location_block(fp, &m, s, i, tttp, ttts, gh, eikonal, dx, dz);
        fclose(fp);
        fprintf(stderr, "wrote %d rows (picks) across %d events to %s\n", total, ne, outf);
    }
    return 0;
}

// ===========================================================================
// EXTENSION NOTES toward global 1c (all events + shared velocity + station cols)
// ---------------------------------------------------------------------------
// Parameter vector ordering (suggested):
//   [ (x,y,z,ot)_0 ... (x,y,z,ot)_{ne-1} | vp_0..vp_{K-1} | vpvs_0..vpvs_{K-1}
//     | pres_0..pres_{S-1} | sres_0..sres_{S-1} ]
// Row set: every P and S pick of every event (global).
//
// - Location columns: as here (central FD via pred_time_pick), placed in the
//   block-diagonal position for the owning event; zero for other events.
// - Station-correction columns: ANALYTIC. dTp/dpres[s] = +1 for a P pick at
//   station s (0 otherwise); dTs/dsres[s] = +1 for an S pick at station s.
//   (Watch the zero-mean / reference-station constraint the MCMC imposes: with
//   scor_flag==0 the corrections are constrained sum=0, so one column is
//   dependent and G^T G is singular unless you drop a column or add the
//   constraint as a pseudo-row.)
// - Velocity columns (vp_k, vpvs_k): perturb m.vp[k] (or m.vpvs[k]) by delta,
//   REBUILD the tables with setup_table_new for the affected phase(s):
//     vp_k   changes both P and S tables (S vel = vp/vpvs) -> rebuild both;
//     vpvs_k changes only the S table -> rebuild S only.
//   Then recompute t_pred for ALL picks and difference. This is the expensive
//   part: 2 table rebuilds per vp column (central diff) x K columns. Cache the
//   base tables and restore by swap, mirroring the mcmc_eq pointer-swap trick.
//   Because velocity enters through find_in_cell's nearest-cell (Voronoi) map,
//   a change to vp[k] affects a contiguous depth range; a full recompute is the
//   safe choice for correctness.
//
// The covariance then follows from the assembled global G:
//     Cm = (G^T Cd^-1 G)^-1 ,  Cd^-1 = diag(1/sigma_row^2)
// with sigma_row taken from the class-appropriate noise reported per row here.
// ===========================================================================
