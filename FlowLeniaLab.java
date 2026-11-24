/*
 * FlowLeniaLab.java
 *
 * A single-file Java Swing laboratory for Flow Lenia, the mass-conserving
 * extension of Lenia introduced by Erwan Plantec, Gautier Hamon, Mayalen
 * Etcheverry, Pierre-Yves Oudeyer, Clement Moulin-Frier and Bert Wang-Chak Chan:
 *
 *   "Flow-Lenia: Towards open-ended evolution in cellular automata through
 *    mass conservation and parameter localization", ALIFE 2023.
 *    https://arxiv.org/abs/2212.07906
 *
 * Lenia (by Bert Wang-Chak Chan) computes a growth field from ring-shaped
 * kernels. Flow Lenia turns that field into a FORCE: matter flows up the
 * growth gradient and away from crowding, and is moved (never created or
 * destroyed) by "reintegration tracking". In the ecosystem mode each cell also
 * carries its own copy of the rule's kernel weights, which travel with the
 * matter, so different species can share one world, and small random
 * mutations let them evolve.
 *
 * This is an independent Java implementation written from the equations in
 * the paper and the authors' published model; it is not a copy of their code.
 *
 * Run:   java FlowLeniaLab.java        (Java 17 or newer)
 *
 * Architect: Mario Gianota
 * Implementation: Claude (Anthropic)
 * Licence: MIT
 */

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.Path2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.*;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.IntStream;

public class FlowLeniaLab extends JFrame {

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            try { UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName()); }
            catch (Exception ignored) { }
            new FlowLeniaLab().setVisible(true);
        });
    }

    // =====================================================================
    //  Model (FFT, rule, world)
    // =====================================================================

    static final class FFT2D {
        final int n;
        private final double[] cos, sin;
        private final int[] rev;
        private final ThreadLocal<double[][]> scratch;
        FFT2D(int n) {
            this.n = n;
            cos = new double[n / 2]; sin = new double[n / 2];
            for (int k = 0; k < n / 2; k++) { cos[k] = Math.cos(2 * Math.PI * k / n); sin[k] = Math.sin(2 * Math.PI * k / n); }
            rev = new int[n];
            int bits = Integer.numberOfTrailingZeros(n);
            for (int i = 0; i < n; i++) rev[i] = Integer.reverse(i) >>> (32 - bits);
            scratch = ThreadLocal.withInitial(() -> new double[2][n]);
        }
        private void fft1d(double[] re, double[] im, boolean inverse) {
            for (int i = 0; i < n; i++) { int j = rev[i]; if (i < j) { double t = re[i]; re[i] = re[j]; re[j] = t; t = im[i]; im[i] = im[j]; im[j] = t; } }
            for (int len = 2; len <= n; len <<= 1) {
                int half = len >> 1, step = n / len;
                for (int i = 0; i < n; i += len)
                    for (int k = 0; k < half; k++) {
                        double wr = cos[k * step], wi = inverse ? sin[k * step] : -sin[k * step];
                        int a = i + k, b = a + half;
                        double xr = re[b] * wr - im[b] * wi, xi = re[b] * wi + im[b] * wr;
                        re[b] = re[a] - xr; im[b] = im[a] - xi; re[a] += xr; im[a] += xi;
                    }
            }
        }
        void forward(double[] re, double[] im) { transform(re, im, false); }
        void inverse(double[] re, double[] im) {
            transform(re, im, true);
            double s = 1.0 / ((double) n * n);
            for (int i = 0; i < re.length; i++) { re[i] *= s; im[i] *= s; }
        }
        private void transform(double[] re, double[] im, boolean inv) {
            IntStream.range(0, n).parallel().forEach(y -> {
                double[][] t = scratch.get();
                System.arraycopy(re, y * n, t[0], 0, n); System.arraycopy(im, y * n, t[1], 0, n);
                fft1d(t[0], t[1], inv);
                System.arraycopy(t[0], 0, re, y * n, n); System.arraycopy(t[1], 0, im, y * n, n);
            });
            IntStream.range(0, n).parallel().forEach(x -> {
                double[][] t = scratch.get();
                for (int y = 0; y < n; y++) { t[0][y] = re[y * n + x]; t[1][y] = im[y * n + x]; }
                fft1d(t[0], t[1], inv);
                for (int y = 0; y < n; y++) { re[y * n + x] = t[0][y]; im[y * n + x] = t[1][y]; }
            });
        }
    }

    /** A Flow Lenia rule: shared kernel shapes plus per-kernel growth parameters. */
    static final class Rule {
        final long seed;
        final int C = 3;
        final int[] c0, c1;          // kernel k reads channel c0[k] and writes to channel c1[k]
        final int K;
        double R;
        final double[] r, m, s, h;
        final double[][] a, b, w;

        /** Random rule, with the same ranges as the reference implementation. */
        Rule(long seed) {
            this.seed = seed;
            int[][] M = {{3, 1, 0}, {0, 3, 1}, {1, 0, 3}};
            List<Integer> src = new ArrayList<>(), dst = new ArrayList<>();
            for (int si = 0; si < 3; si++) for (int t = 0; t < 3; t++) for (int q = 0; q < M[si][t]; q++) { src.add(si); dst.add(t); }
            K = src.size();
            c0 = src.stream().mapToInt(Integer::intValue).toArray();
            c1 = dst.stream().mapToInt(Integer::intValue).toArray();
            Random rng = new Random(mix(seed));
            R = uni(rng, 2, 25);
            r = new double[K]; m = new double[K]; s = new double[K]; h = new double[K];
            a = new double[K][3]; b = new double[K][3]; w = new double[K][3];
            for (int k = 0; k < K; k++) {
                r[k] = uni(rng, 0.2, 1.0);
                m[k] = uni(rng, 0.05, 0.5);
                s[k] = uni(rng, 0.001, 0.18);
                h[k] = uni(rng, 0.01, 1.0);
                for (int j = 0; j < 3; j++) { a[k][j] = uni(rng, 0, 1); b[k][j] = uni(rng, 0.001, 1); w[k][j] = uni(rng, 0.01, 0.5); }
            }
        }
        /** SplitMix64 scramble, so neighbouring seeds give unrelated rules. */
        static long mix(long z) { z += 0x9e3779b97f4a7c15L; z = (z ^ (z >>> 30)) * 0xbf58476d1ce4e5b9L; z = (z ^ (z >>> 27)) * 0x94d049bb133111ebL; return z ^ (z >>> 31); }

        static double uni(Random r, double lo, double hi) { return lo + (hi - lo) * r.nextDouble(); }

        static double sigmoid(double x) { return 0.5 * (Math.tanh(x / 2) + 1); }

        /** Kernel k's radial profile at distance d (cells). */
        double kernel(int k, double d) {
            double D = d / ((R + 15) * r[k]);
            double v = 0;
            for (int j = 0; j < 3; j++) v += b[k][j] * Math.exp(-(D - a[k][j]) * (D - a[k][j]) / w[k][j]);
            return sigmoid(-(D - 1) * 10) * v;
        }
    }

    static final class World {
        final int n, C = 3;
        final FFT2D fft;
        Rule rule;
        double[][] A;                // [channel][cell]
        double[][] P;                // [kernel][cell], or null in single-rule mode
        double dt = 0.2, sigma = 0.65;
        final int dd = 5;
        double theta;                // alpha = (A / theta)^2
        boolean torus = true;
        long gen = 0;
        final double[][] kre, kim;
        final double[] re, im;
        final double[][] U;          // growth field per channel
        double[][] nA;
        double[][] nP;
        final double[] tot;
        final int[] chosen;
        final Random rng = new Random(1);
        double mutationProb = 0, mutationSize = 20, mutationStrength = 1.0;

        World(int n, Rule rule, boolean multi) {
            this.n = n;
            fft = new FFT2D(n);
            int nn = n * n;
            A = new double[C][nn]; nA = new double[C][nn];
            U = new double[C][nn];
            re = new double[nn]; im = new double[nn];
            tot = new double[nn]; chosen = new int[nn];
            kre = new double[rule.K][nn]; kim = new double[rule.K][nn];
            setRule(rule, multi);
        }

        void setRule(Rule rule, boolean multi) {
            this.rule = rule;
            int nn = n * n;
            theta = multi ? 2.0 : C;
            P = multi ? new double[rule.K][nn] : null;
            nP = multi ? new double[rule.K][nn] : null;
            for (int k = 0; k < rule.K; k++) {
                double sum = 0;
                for (int y = 0; y < n; y++) {
                    int dy = y < n / 2 ? y : y - n;
                    for (int x = 0; x < n; x++) {
                        int dx = x < n / 2 ? x : x - n;
                        double v = rule.kernel(k, Math.sqrt(dx * dx + dy * dy));
                        kre[k][y * n + x] = v; sum += v;
                    }
                }
                for (int i = 0; i < nn; i++) { kre[k][i] /= sum; kim[k][i] = 0; }
                fft.forward(kre[k], kim[k]);
            }
        }

        static double growth(double u, double m, double s) { double z = (u - m) / s; return Math.exp(-z * z / 2) * 2 - 1; }

        void step() {
            int nn = n * n;
            Rule R = rule;
            for (double[] u : U) Arrays.fill(u, 0);
            // Forward FFT of each channel
            double[][] fre = new double[C][], fim = new double[C][];
            for (int c = 0; c < C; c++) {
                fre[c] = A[c].clone(); fim[c] = new double[nn];
                fft.forward(fre[c], fim[c]);
            }
            // Each kernel: convolve, apply growth, weight by h (or by embedded P), add to target channel
            for (int k = 0; k < R.K; k++) {
                double[] ar = fre[R.c0[k]], ai = fim[R.c0[k]], kr = kre[k], ki = kim[k];
                for (int i = 0; i < nn; i++) { re[i] = ar[i] * kr[i] - ai[i] * ki[i]; im[i] = ar[i] * ki[i] + ai[i] * kr[i]; }
                fft.inverse(re, im);
                double m = R.m[k], s = R.s[k], h = R.h[k];
                double[] out = U[R.c1[k]];
                double[] pk = P == null ? null : P[k];
                for (int i = 0; i < nn; i++) out[i] += growth(re[i], m, s) * (pk == null ? h : pk[i]);
            }
            // Flow: F = grad(U) (1 - alpha) - grad(sum A) alpha, then reintegration tracking
            double[] Asum = new double[nn];
            for (int c = 0; c < C; c++) for (int i = 0; i < nn; i++) Asum[i] += A[c][i];
            double[][] gA = sobel(Asum);
            for (double[] x : nA) Arrays.fill(x, 0);
            if (P != null) { Arrays.fill(tot, 0); Arrays.fill(chosen, -1); }
            double ma = dd - sigma;
            double inv4s2 = 1.0 / (4 * sigma * sigma), cap = Math.min(1, 2 * sigma);
            double[][] gU = new double[C][];
            double[][][] gUc = new double[C][][];
            for (int c = 0; c < C; c++) gUc[c] = sobel(U[c]);
            // Per-source accumulation buffer for P mixing
            int[] dstIdx = new int[64]; double[] dstW = new double[64];
            for (int i = 0; i < nn; i++) {
                int y = i / n, x = i % n;
                int cnt = 0;
                for (int c = 0; c < C; c++) {
                    double a = A[c][i];
                    if (a == 0) continue;
                    double alpha = a / theta; alpha = Math.min(1, alpha * alpha);
                    double fy = gUc[c][0][i] * (1 - alpha) - gA[0][i] * alpha;   // along axis 0 (rows)
                    double fx = gUc[c][1][i] * (1 - alpha) - gA[1][i] * alpha;   // along axis 1 (columns)
                    double muy = y + 0.5 + clamp(dt * fy, -ma, ma);
                    double mux = x + 0.5 + clamp(dt * fx, -ma, ma);
                    if (!torus) { muy = clamp(muy, sigma, n - sigma); mux = clamp(mux, sigma, n - sigma); }
                    int y0 = (int) Math.floor(muy - sigma), y1 = (int) Math.floor(muy + sigma);
                    int x0 = (int) Math.floor(mux - sigma), x1 = (int) Math.floor(mux + sigma);
                    for (int py = y0; py <= y1; py++) {
                        double oy = Math.min(py + 1, muy + sigma) - Math.max(py, muy - sigma);
                        if (oy <= 0) continue;
                        oy = Math.min(oy, cap);
                        for (int px = x0; px <= x1; px++) {
                            double ox = Math.min(px + 1, mux + sigma) - Math.max(px, mux - sigma);
                            if (ox <= 0) continue;
                            ox = Math.min(ox, cap);
                            int d = Math.floorMod(py, n) * n + Math.floorMod(px, n);
                            double moved = a * oy * ox * inv4s2;
                            nA[c][d] += moved;
                            if (P != null) {
                                int q = 0;
                                while (q < cnt && dstIdx[q] != d) q++;
                                if (q == cnt) { dstIdx[cnt] = d; dstW[cnt] = 0; cnt++; }
                                dstW[q] += moved;
                            }
                        }
                    }
                }
                if (P != null)
                    for (int q = 0; q < cnt; q++) {
                        int d = dstIdx[q];
                        tot[d] += dstW[q];
                        if (rng.nextDouble() * tot[d] < dstW[q]) chosen[d] = i;   // mass-weighted random pick
                    }
            }
            double[][] t = A; A = nA; nA = t;
            if (P != null) {
                for (int k = 0; k < rule.K; k++) {
                    double[] src = P[k], dst = nP[k];
                    for (int d = 0; d < nn; d++) dst[d] = chosen[d] >= 0 ? src[chosen[d]] : 0;
                }
                double[][] tp = P; P = nP; nP = tp;
                if (mutationProb > 0 && rng.nextDouble() < mutationProb) mutate();
            }
            gen++;
        }

        /** Add the same random change to the embedded parameters of a square region (a "beam" mutation). */
        void mutate() {
            int sz = (int) mutationSize;
            int cy = rng.nextInt(n), cx = rng.nextInt(n);
            double[] dv = new double[rule.K];
            for (int k = 0; k < dv.length; k++) dv[k] = rng.nextGaussian() * mutationStrength;
            for (int y = 0; y < sz; y++) for (int x = 0; x < sz; x++) {
                int i = Math.floorMod(cy + y, n) * n + Math.floorMod(cx + x, n);
                for (int k = 0; k < dv.length; k++) P[k][i] += dv[k];
            }
        }

        static double clamp(double v, double lo, double hi) { return v < lo ? lo : v > hi ? hi : v; }

        /** Sobel gradients {d/drow, d/dcol}, matching the reference orientation, with wrap-around. */
        double[][] sobel(double[] f) {
            int nn = n * n;
            double[] gy = new double[nn], gx = new double[nn];
            IntStream.range(0, n).parallel().forEach(y -> {
                int ym = (y - 1 + n) % n, yp = (y + 1) % n;
                for (int x = 0; x < n; x++) {
                    int xm = (x - 1 + n) % n, xp = (x + 1) % n;
                    double a = f[ym * n + xm], b = f[ym * n + x], c = f[ym * n + xp];
                    double d = f[y * n + xm], e = f[y * n + xp];
                    double g = f[yp * n + xm], hh = f[yp * n + x], ii = f[yp * n + xp];
                    gx[y * n + x] = (c + 2 * e + ii) - (a + 2 * d + g);
                    gy[y * n + x] = (g + 2 * hh + ii) - (a + 2 * b + c);
                }
            });
            return new double[][]{gy, gx};
        }

        double mass() { double s = 0; for (double[] c : A) for (double v : c) s += v; return s; }

        void seedBlock(int cy, int cx, int size, Random r, double[] pvec) {
            for (int y = 0; y < size; y++) for (int x = 0; x < size; x++) {
                int i = Math.floorMod(cy - size / 2 + y, n) * n + Math.floorMod(cx - size / 2 + x, n);
                for (int c = 0; c < C; c++) A[c][i] = r.nextDouble();
                if (P != null && pvec != null) for (int k = 0; k < rule.K; k++) P[k][i] = pvec[k];
            }
        }
    }


    // =====================================================================
    //  Species colours
    // =====================================================================

    static final int HUE_BINS = 24;

    /** Fixed random projection of a cell's embedded parameters onto a colour wheel. */
    static final class Palette {
        final double[] px, py;
        Palette(int K) {
            Random r = new Random(12345);
            px = new double[K]; py = new double[K];
            for (int k = 0; k < K; k++) { px[k] = r.nextGaussian(); py[k] = r.nextGaussian(); }
        }
        double hue(double[][] P, int i) {
            double x = 0, y = 0;
            for (int k = 0; k < px.length; k++) { double v = P[k][i] - 0.5; x += px[k] * v; y += py[k] * v; }
            return (Math.atan2(y, x) / (2 * Math.PI) + 1) % 1.0;
        }
        double hue(double[] p) {
            double x = 0, y = 0;
            for (int k = 0; k < px.length; k++) { double v = p[k] - 0.5; x += px[k] * v; y += py[k] * v; }
            return (Math.atan2(y, x) / (2 * Math.PI) + 1) % 1.0;
        }
    }

    enum View {
        SPECIES("Species"), MATTER("Matter (3 channels)");
        final String label;
        View(String l) { label = l; }
        @Override public String toString() { return label; }
    }

    enum Mode {
        ECOSYSTEM("Ecosystem (many species)"), SINGLE("Single rule");
        final String label;
        Mode(String l) { label = l; }
        @Override public String toString() { return label; }
    }

    // =====================================================================
    //  Application state
    // =====================================================================

    private volatile World world;
    private volatile Palette palette;
    private final ConcurrentLinkedQueue<Runnable> tasks = new ConcurrentLinkedQueue<>();
    private volatile boolean running = true, singleStep = false;
    private volatile int stepsPerFrame = 2;
    private volatile View view = View.SPECIES;
    private volatile Mode mode = Mode.ECOSYSTEM;
    private volatile BufferedImage frame;
    private volatile double stepsPerSecond, lastMass, diversity;
    private volatile int creatures;
    private volatile int speciesCount = 6;
    private volatile boolean evolution = false;
    private volatile double mutationRate = 0.02, mutationStrength = 0.15;
    private final Random rng = new Random();

    /** History of each colour's share of the world's matter, for the history strip. */
    private final ArrayList<double[]> history = new ArrayList<>();
    private final ArrayList<Double> diversityHistory = new ArrayList<>();
    private final Object histLock = new Object();

    private final WorldView worldView = new WorldView();
    private final HistoryStrip historyStrip = new HistoryStrip();
    private final KernelPlot kernelPlot = new KernelPlot();
    private final JLabel status = new JLabel(" ");
    private final JLabel ruleInfo = new JLabel(" ");

    private final JComboBox<Mode> modeBox = new JComboBox<>(Mode.values());
    private final JComboBox<View> viewBox = new JComboBox<>(View.values());
    private final JComboBox<Integer> sizeBox = new JComboBox<>(new Integer[]{128, 256});
    private final JTextField seedField = new JTextField("5", 8);
    private final ParamSlider speciesSlider = new ParamSlider("Starting species", 1, 12, 6, "%.0f");
    private final ParamSlider dtSlider = new ParamSlider("Time step dt", 0.05, 0.5, 0.2, "%.2f");
    private final ParamSlider rateSlider = new ParamSlider("Mutation rate (per step)", 0, 0.1, 0.02, "%.3f");
    private final ParamSlider strengthSlider = new ParamSlider("Mutation strength", 0.01, 1.0, 0.15, "%.2f");
    private final ParamSlider speedSlider = new ParamSlider("Steps per frame", 1, 20, 2, "%.0f");
    private final JCheckBox evoBox = new JCheckBox("Evolution: random mutations");
    private final JToggleButton playButton = new JToggleButton("Pause", true);

    public FlowLeniaLab() {
        super("Flow Lenia Lab — an ecosystem of artificial creatures");
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        Rule r = new Rule(5);
        world = new World(128, r, true);
        palette = new Palette(r.K);
        seedWorld(world);
        buildUI();
        updateRuleInfo(r);

        Thread worker = new Thread(this::loop, "flow-lenia-sim");
        worker.setDaemon(true);
        worker.start();
        new javax.swing.Timer(60, e -> { worldView.repaint(); historyStrip.repaint(); updateStatus(); }).start();
        pack();
        setLocationRelativeTo(null);
    }

    // =====================================================================
    //  Simulation loop
    // =====================================================================

    private void loop() {
        World current = null;
        int[] pixels = null;
        BufferedImage[] bufs = new BufferedImage[2];
        int bi = 0, frameNo = 0;
        long lastT = System.nanoTime(), stepsSince = 0;
        while (true) {
            long start = System.nanoTime();
            Runnable t;
            while ((t = tasks.poll()) != null) t.run();
            World w = world;
            if (w != current) {
                current = w;
                bufs[0] = new BufferedImage(w.n, w.n, BufferedImage.TYPE_INT_RGB);
                bufs[1] = new BufferedImage(w.n, w.n, BufferedImage.TYPE_INT_RGB);
                pixels = new int[w.n * w.n];
            }
            w.mutationProb = evolution && w.P != null ? mutationRate : 0;
            w.mutationStrength = mutationStrength;
            boolean advance = running || singleStep;
            if (advance) {
                int steps = singleStep && !running ? 1 : stepsPerFrame;
                singleStep = false;
                for (int k = 0; k < steps; k++) {
                    w.step();
                    stepsSince++;
                    if (w.gen % 10 == 0) recordHistory(w);
                }
            }
            render(w, pixels);
            BufferedImage img = bufs[bi];
            bi ^= 1;
            img.setRGB(0, 0, w.n, w.n, pixels, 0, w.n);
            frame = img;
            if (frameNo % 10 == 0) { lastMass = w.mass(); creatures = countCreatures(w); }
            frameNo++;

            long now = System.nanoTime();
            if (now - lastT > 500_000_000L) { stepsPerSecond = stepsSince * 1e9 / (now - lastT); stepsSince = 0; lastT = now; }
            long sleep = 25 - (System.nanoTime() - start) / 1_000_000;
            try { Thread.sleep(Math.max(2, sleep)); } catch (InterruptedException ex) { return; }
        }
    }

    private void render(World w, int[] px) {
        int nn = w.n * w.n;
        if (view == View.SPECIES && w.P != null) {
            Palette pal = palette;
            for (int i = 0; i < nn; i++) {
                double m = w.A[0][i] + w.A[1][i] + w.A[2][i];
                if (m < 0.004) { px[i] = 0x05050a; continue; }
                float br = (float) (1 - Math.exp(-1.3 * m));
                px[i] = Color.HSBtoRGB((float) pal.hue(w.P, i), 0.78f, br);
            }
        } else {
            for (int i = 0; i < nn; i++) {
                int r = (int) (255 * (1 - Math.exp(-1.5 * w.A[0][i])));
                int g = (int) (255 * (1 - Math.exp(-1.5 * w.A[1][i])));
                int b = (int) (255 * (1 - Math.exp(-1.5 * w.A[2][i])));
                px[i] = (r << 16) | (g << 8) | b;
            }
        }
    }

    /** Share of the world's matter in each colour band, plus an effective number of species. */
    private void recordHistory(World w) {
        int nn = w.n * w.n;
        double[] bins = new double[HUE_BINS];
        double tot = 0;
        if (w.P != null) {
            Palette pal = palette;
            for (int i = 0; i < nn; i += 2) {
                double m = w.A[0][i] + w.A[1][i] + w.A[2][i];
                if (m < 0.01) continue;
                int b = Math.min(HUE_BINS - 1, (int) (pal.hue(w.P, i) * HUE_BINS));
                bins[b] += m;
                tot += m;
            }
        } else {
            for (int c = 0; c < 3; c++) { double s = 0; for (int i = 0; i < nn; i += 2) s += w.A[c][i]; bins[c * 8] = s; tot += s; }
        }
        double H = 0;
        if (tot > 0) for (int b = 0; b < HUE_BINS; b++) { bins[b] /= tot; if (bins[b] > 0) H -= bins[b] * Math.log(bins[b]); }
        double div = Math.exp(H);
        diversity = div;
        synchronized (histLock) {
            history.add(bins);
            diversityHistory.add(div);
            if (history.size() > 2000) {
                ArrayList<double[]> h2 = new ArrayList<>();
                ArrayList<Double> d2 = new ArrayList<>();
                for (int i = 0; i < history.size(); i += 2) { h2.add(history.get(i)); d2.add(diversityHistory.get(i)); }
                history.clear(); history.addAll(h2);
                diversityHistory.clear(); diversityHistory.addAll(d2);
            }
        }
    }

    /** Count separate bodies of matter (connected regions holding a meaningful amount of mass). */
    private static int countCreatures(World w) {
        int n = w.n, nn = n * n;
        boolean[] seen = new boolean[nn];
        int[] stack = new int[nn];
        int count = 0;
        for (int s = 0; s < nn; s++) {
            if (seen[s] || w.A[0][s] + w.A[1][s] + w.A[2][s] < 0.05) continue;
            int sp = 0; stack[sp++] = s; seen[s] = true;
            double mass = 0;
            while (sp > 0) {
                int c = stack[--sp];
                mass += w.A[0][c] + w.A[1][c] + w.A[2][c];
                int x = c % n, y = c / n;
                int[] nb = {y * n + (x + 1) % n, y * n + (x - 1 + n) % n, ((y + 1) % n) * n + x, ((y - 1 + n) % n) * n + x};
                for (int j : nb) if (!seen[j] && w.A[0][j] + w.A[1][j] + w.A[2][j] >= 0.05) { seen[j] = true; stack[sp++] = j; }
            }
            if (mass > 5) count++;
        }
        return count;
    }

    // =====================================================================
    //  Actions
    // =====================================================================

    private void seedWorld(World w) {
        Random r = new Random(rng.nextLong());
        int size = Math.max(12, w.n / 6);
        if (w.P == null) {
            w.seedBlock(w.n / 2, w.n / 2, Math.max(24, w.n * 40 / 128), r, null);
        } else {
            int sp = speciesCount;
            for (int s = 0; s < sp; s++) {
                double[] p = new double[w.rule.K];
                for (int k = 0; k < p.length; k++) p[k] = r.nextDouble();
                w.seedBlock(r.nextInt(w.n), r.nextInt(w.n), size, r, p);
            }
        }
    }

    private void restart(long seed) {
        int n = (Integer) sizeBox.getSelectedItem();
        boolean multi = mode == Mode.ECOSYSTEM;
        double dt = dtSlider.value();
        tasks.add(() -> {
            Rule r = new Rule(seed);
            World w = new World(n, r, multi);
            w.dt = dt;
            palette = new Palette(r.K);
            seedWorld(w);
            world = w;
            synchronized (histLock) { history.clear(); diversityHistory.clear(); }
            SwingUtilities.invokeLater(() -> { updateRuleInfo(r); kernelPlot.repaint(); });
        });
    }

    private long readSeed() {
        try { return Long.parseLong(seedField.getText().trim()); }
        catch (NumberFormatException e) { seedField.setText("5"); return 5; }
    }

    private void updateRuleInfo(Rule r) {
        ruleInfo.setText(String.format(Locale.ROOT, "<html><div style='width:250px'>Rule %d: %d kernels linking 3 channels, base radius R = %.1f</div></html>", r.seed, r.K, r.R));
    }

    private void updateStatus() {
        World w = world;
        String div = w.P != null ? String.format(Locale.ROOT, "   diversity %.1f", diversity) : "";
        status.setText(String.format(Locale.ROOT, "step %,8d   grid %d×%d   %5.1f steps/s   total mass %.4f   bodies %d%s%s",
                w.gen, w.n, w.n, stepsPerSecond, lastMass, creatures, div, evolution && w.P != null ? "   evolving" : ""));
    }

    private void savePng() {
        BufferedImage img = frame;
        if (img == null) return;
        JFileChooser fc = new JFileChooser();
        fc.setSelectedFile(new File("flow-lenia.png"));
        if (fc.showSaveDialog(this) != JFileChooser.APPROVE_OPTION) return;
        File f = fc.getSelectedFile();
        if (!f.getName().toLowerCase().endsWith(".png")) f = new File(f.getPath() + ".png");
        BufferedImage out = new BufferedImage(1024, 1024, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = out.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
        g.drawImage(img, 0, 0, 1024, 1024, null);
        g.dispose();
        try { ImageIO.write(out, "png", f); }
        catch (Exception ex) { JOptionPane.showMessageDialog(this, "Could not save: " + ex.getMessage()); }
    }

    private void showAbout() {
        String html = """
            <html><body style='width:500px; font-family:sans-serif'>
            <h2>Flow Lenia</h2>
            <p><b>Lenia</b>, created by Bert Wang-Chak Chan, grows smooth artificial creatures from
            a simple rule: each cell senses a soft ring of neighbours and grows or fades depending
            on what it senses. Its creatures are beautiful, but each world has a single rule, so
            different species can never meet, and life can appear from nothing or vanish.</p>
            <p><b>Flow Lenia</b> (Plantec, Hamon, Etcheverry, Oudeyer, Moulin-Frier and Chan, 2023)
            changes two things:</p>
            <ul>
            <li><b>Matter is conserved.</b> Lenia's growth signal becomes a force. Matter flows towards
            better growing conditions and away from crowding, and is moved, never created or
            destroyed. Watch the total mass in the status bar: it never changes.</li>
            <li><b>The rule travels with the matter.</b> In the ecosystem mode, every cell carries
            its own weights for the rule's kernels. When matter flows, its weights go with it. So
            each colour you see is a different species with its own rule, and species can meet,
            push, merge and compete in one world.</li>
            </ul>
            <p>Switch on <b>Evolution</b> and every so often a random patch of the world has its
            weights nudged. Most mutants vanish; some spread. The <b>history strip</b> records each
            colour's share of the world over time, so you can watch lineages rise and fall.</p>
            <p><b>Matter view</b> shows the three channels of matter as red, green and blue.
            <b>New rule</b> draws a fresh random set of kernels, as in the paper's random search:
            almost every rule produces creatures, but each looks different.</p>
            <p>Based on <i>Flow-Lenia: Towards open-ended evolution in cellular automata through
            mass conservation and parameter localization</i> (ALIFE 2023, arXiv:2212.07906), and on
            Lenia by Bert Wang-Chak Chan. An independent implementation, written from the paper.</p>
            <p style='color:#666'>Architect: Mario Gianota · Implementation: Claude (Anthropic) · MIT licence</p>
            </body></html>""";
        JOptionPane.showMessageDialog(this, new JLabel(html), "About Flow Lenia Lab", JOptionPane.PLAIN_MESSAGE);
    }

    // =====================================================================
    //  UI
    // =====================================================================

    private void buildUI() {
        JPanel controls = new JPanel();
        controls.setLayout(new BoxLayout(controls, BoxLayout.Y_AXIS));
        controls.setBorder(new EmptyBorder(8, 8, 8, 14));

        JPanel worldP = section("World");
        worldP.add(row("Mode: ", modeBox));
        modeBox.addActionListener(e -> {
            mode = (Mode) modeBox.getSelectedItem();
            if (mode == Mode.SINGLE) viewBox.setSelectedItem(View.MATTER);
            else viewBox.setSelectedItem(View.SPECIES);
            restart(readSeed());
        });
        JPanel seedRow = row("Rule seed: ", seedField);
        JButton useSeed = new JButton("Start");
        useSeed.addActionListener(e -> restart(readSeed()));
        seedRow.add(Box.createHorizontalStrut(6));
        seedRow.add(useSeed);
        worldP.add(seedRow);
        ruleInfo.setFont(ruleInfo.getFont().deriveFont(11f));
        ruleInfo.setAlignmentX(LEFT_ALIGNMENT);
        worldP.add(ruleInfo);
        kernelPlot.setAlignmentX(LEFT_ALIGNMENT);
        worldP.add(Box.createVerticalStrut(4));
        worldP.add(kernelPlot);
        speciesSlider.onChange(() -> speciesCount = (int) Math.round(speciesSlider.value()));
        worldP.add(speciesSlider);
        JPanel b1 = new JPanel(new GridLayout(1, 2, 4, 4));
        JButton newRule = new JButton("New rule");
        JButton reseed = new JButton("Reseed");
        newRule.setToolTipText("Draw a random new rule, as in the paper's random search");
        reseed.setToolTipText("Keep the rule, scatter fresh species");
        newRule.addActionListener(e -> { long s = Math.abs(rng.nextInt(100000)); seedField.setText(String.valueOf(s)); restart(s); });
        reseed.addActionListener(e -> restart(readSeed()));
        b1.add(newRule); b1.add(reseed);
        b1.setAlignmentX(LEFT_ALIGNMENT);
        worldP.add(b1);
        JPanel opts = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 4));
        opts.add(new JLabel("Grid: "));
        sizeBox.addActionListener(e -> restart(readSeed()));
        opts.add(sizeBox);
        opts.add(Box.createHorizontalStrut(10));
        opts.add(new JLabel("Show: "));
        viewBox.addActionListener(e -> view = (View) viewBox.getSelectedItem());
        opts.add(viewBox);
        opts.setAlignmentX(LEFT_ALIGNMENT);
        worldP.add(opts);
        controls.add(worldP);

        JPanel evo = section("Evolution");
        evoBox.setAlignmentX(LEFT_ALIGNMENT);
        evoBox.addActionListener(e -> evolution = evoBox.isSelected());
        evo.add(evoBox);
        rateSlider.onChange(() -> mutationRate = rateSlider.value());
        strengthSlider.onChange(() -> mutationStrength = strengthSlider.value());
        evo.add(rateSlider);
        evo.add(strengthSlider);
        JLabel hl = new JLabel("History: each colour's share over time");
        hl.setFont(hl.getFont().deriveFont(11f));
        hl.setAlignmentX(LEFT_ALIGNMENT);
        evo.add(hl);
        historyStrip.setAlignmentX(LEFT_ALIGNMENT);
        evo.add(historyStrip);
        controls.add(evo);

        JPanel sim = section("Simulation");
        dtSlider.onChange(() -> { double v = dtSlider.value(); tasks.add(() -> world.dt = v); });
        speedSlider.onChange(() -> stepsPerFrame = (int) Math.round(speedSlider.value()));
        sim.add(dtSlider);
        sim.add(speedSlider);
        JPanel b2 = new JPanel(new GridLayout(2, 2, 4, 4));
        JButton stepBtn = new JButton("Step");
        JButton saveBtn = new JButton("Save PNG");
        JButton aboutBtn = new JButton("About");
        playButton.addActionListener(e -> { running = playButton.isSelected(); playButton.setText(running ? "Pause" : "Play"); });
        stepBtn.addActionListener(e -> singleStep = true);
        saveBtn.addActionListener(e -> savePng());
        aboutBtn.addActionListener(e -> showAbout());
        b2.add(playButton); b2.add(stepBtn); b2.add(saveBtn); b2.add(aboutBtn);
        b2.setAlignmentX(LEFT_ALIGNMENT);
        sim.add(Box.createVerticalStrut(4));
        sim.add(b2);
        controls.add(sim);

        JLabel hint = new JLabel("<html><i>Left-click: drop in a new species (or matter)<br>Right-drag: remove matter</i></html>");
        hint.setBorder(new EmptyBorder(6, 4, 0, 0));
        hint.setAlignmentX(LEFT_ALIGNMENT);
        controls.add(hint);
        controls.add(Box.createVerticalGlue());

        JScrollPane scroll = new JScrollPane(controls, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
                ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
        scroll.setBorder(null);
        scroll.setPreferredSize(new Dimension(400, 820));
        scroll.getVerticalScrollBar().setUnitIncrement(16);

        status.setBorder(new EmptyBorder(4, 10, 4, 10));
        status.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));

        getContentPane().setLayout(new BorderLayout());
        getContentPane().add(scroll, BorderLayout.WEST);
        getContentPane().add(worldView, BorderLayout.CENTER);
        getContentPane().add(status, BorderLayout.SOUTH);
    }

    private static JPanel row(String label, JComponent c) {
        JPanel p = new JPanel(new FlowLayout(FlowLayout.LEFT, 0, 2));
        p.add(new JLabel(label));
        p.add(c);
        p.setAlignmentX(LEFT_ALIGNMENT);
        return p;
    }

    private static JPanel section(String title) {
        JPanel p = new JPanel();
        p.setLayout(new BoxLayout(p, BoxLayout.Y_AXIS));
        p.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), title, TitledBorder.LEFT, TitledBorder.TOP),
                new EmptyBorder(4, 6, 6, 6)));
        p.setAlignmentX(LEFT_ALIGNMENT);
        return p;
    }

    // =====================================================================
    //  Widgets
    // =====================================================================

    static final class ParamSlider extends JPanel {
        private static final int TICKS = 10000;
        private final double min, max;
        private final String fmt;
        private final JSlider slider = new JSlider(0, TICKS);
        private final JLabel valueLabel = new JLabel();
        private Runnable listener = () -> { };

        ParamSlider(String name, double min, double max, double initial, String fmt) {
            this.min = min; this.max = max; this.fmt = fmt;
            setLayout(new BorderLayout());
            setAlignmentX(LEFT_ALIGNMENT);
            setMaximumSize(new Dimension(Integer.MAX_VALUE, 52));
            JPanel top = new JPanel(new BorderLayout());
            JLabel nameLabel = new JLabel(name);
            nameLabel.setFont(nameLabel.getFont().deriveFont(12f));
            valueLabel.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
            top.add(nameLabel, BorderLayout.WEST);
            top.add(valueLabel, BorderLayout.EAST);
            add(top, BorderLayout.NORTH);
            add(slider, BorderLayout.CENTER);
            slider.addChangeListener(e -> { refresh(); listener.run(); });
            setValue(initial);
        }
        double value() { return min + (max - min) * slider.getValue() / TICKS; }
        void setValue(double v) { slider.setValue((int) Math.round((v - min) / (max - min) * TICKS)); refresh(); }
        void onChange(Runnable r) { listener = r; }
        private void refresh() { valueLabel.setText(String.format(Locale.ROOT, fmt, value())); }
    }

    final class WorldView extends JPanel {
        WorldView() {
            setPreferredSize(new Dimension(800, 800));
            setBackground(new Color(0x08, 0x08, 0x0e));
            MouseAdapter ma = new MouseAdapter() {
                @Override public void mousePressed(MouseEvent e) { handle(e, true); }
                @Override public void mouseDragged(MouseEvent e) { if (SwingUtilities.isRightMouseButton(e)) handle(e, false); }
            };
            addMouseListener(ma);
            addMouseMotionListener(ma);
            setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
        }

        private int side() { return Math.min(getWidth(), getHeight()) - 16; }

        private void handle(MouseEvent e, boolean press) {
            int side = side();
            int ox = (getWidth() - side) / 2, oy = (getHeight() - side) / 2;
            World w = world;
            int gx = (int) ((e.getX() - ox) * (double) w.n / side), gy = (int) ((e.getY() - oy) * (double) w.n / side);
            if (gx < 0 || gy < 0 || gx >= w.n || gy >= w.n) return;
            if (SwingUtilities.isRightMouseButton(e)) {
                int rad = Math.max(3, w.n / 40);
                tasks.add(() -> {
                    if (world != w) return;
                    for (int dy = -rad; dy <= rad; dy++) for (int dx = -rad; dx <= rad; dx++) {
                        if (dx * dx + dy * dy > rad * rad) continue;
                        int i = Math.floorMod(gy + dy, w.n) * w.n + Math.floorMod(gx + dx, w.n);
                        for (int c = 0; c < 3; c++) w.A[c][i] = 0;
                    }
                });
            } else if (press) {
                tasks.add(() -> {
                    if (world != w) return;
                    Random r = new Random(rng.nextLong());
                    double[] p = null;
                    if (w.P != null) { p = new double[w.rule.K]; for (int k = 0; k < p.length; k++) p[k] = r.nextDouble(); }
                    w.seedBlock(gy, gx, Math.max(12, w.n / 6), r, p);
                });
            }
        }

        @Override protected void paintComponent(Graphics g0) {
            super.paintComponent(g0);
            Graphics2D g = (Graphics2D) g0;
            BufferedImage img = frame;
            if (img == null) return;
            int side = side();
            g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);
            g.drawImage(img, (getWidth() - side) / 2, (getHeight() - side) / 2, side, side, null);
        }
    }

    /** Each column is a moment in time; its colours show each species' share of the world's matter. */
    final class HistoryStrip extends JPanel {
        HistoryStrip() {
            setPreferredSize(new Dimension(300, 110));
            setMaximumSize(new Dimension(Integer.MAX_VALUE, 110));
            setBackground(new Color(0x14, 0x14, 0x1c));
        }

        @Override protected void paintComponent(Graphics g0) {
            super.paintComponent(g0);
            Graphics2D g = (Graphics2D) g0;
            List<double[]> h;
            List<Double> d;
            synchronized (histLock) { h = new ArrayList<>(history); d = new ArrayList<>(diversityHistory); }
            int w = getWidth(), ht = getHeight() - 14;
            if (h.isEmpty()) { g.setColor(Color.GRAY); g.drawString("waiting…", 8, ht / 2); return; }
            boolean species = world.P != null;
            int cols = h.size();
            for (int c = 0; c < cols; c++) {
                int x0 = (int) ((double) c / cols * w), x1 = (int) ((double) (c + 1) / cols * w);
                double[] bins = h.get(c);
                double acc = 0;
                for (int b = 0; b < HUE_BINS; b++) {
                    if (bins[b] <= 0) continue;
                    int y0 = (int) (acc * ht), y1 = (int) ((acc + bins[b]) * ht);
                    acc += bins[b];
                    Color col = species ? Color.getHSBColor((b + 0.5f) / HUE_BINS, 0.78f, 0.95f)
                            : (b == 0 ? new Color(230, 70, 70) : b == 8 ? new Color(70, 220, 90) : new Color(80, 120, 240));
                    g.setColor(col);
                    g.fillRect(x0, y0, Math.max(1, x1 - x0), Math.max(1, y1 - y0));
                }
            }
            g.setColor(new Color(210, 210, 210));
            g.setFont(g.getFont().deriveFont(10f));
            g.drawString(species ? String.format(Locale.ROOT, "diversity %.1f  (effective number of colours)", d.get(d.size() - 1)) : "share of each matter channel", 2, getHeight() - 2);
        }
    }

    /** Radial profiles of the rule's kernels, coloured by the channel each one feeds. */
    final class KernelPlot extends JPanel {
        KernelPlot() {
            setPreferredSize(new Dimension(300, 90));
            setMaximumSize(new Dimension(Integer.MAX_VALUE, 90));
            setBackground(new Color(0x14, 0x14, 0x1c));
        }

        @Override protected void paintComponent(Graphics g0) {
            super.paintComponent(g0);
            Graphics2D g = (Graphics2D) g0;
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            Rule r = world.rule;
            int w = getWidth(), h = getHeight() - 14;
            double maxD = (r.R + 15);
            Color[] cols = {new Color(240, 90, 90, 200), new Color(90, 220, 110, 200), new Color(100, 140, 255, 200)};
            double vmax = 1e-9;
            for (int k = 0; k < r.K; k++) for (int i = 0; i <= 100; i++) vmax = Math.max(vmax, r.kernel(k, maxD * i / 100.0));
            g.setStroke(new BasicStroke(1.3f));
            for (int k = 0; k < r.K; k++) {
                Path2D p = new Path2D.Double();
                for (int i = 0; i <= 100; i++) {
                    double x = 4 + (w - 8) * i / 100.0, y = h - 2 - r.kernel(k, maxD * i / 100.0) / vmax * (h - 6);
                    if (i == 0) p.moveTo(x, y); else p.lineTo(x, y);
                }
                g.setColor(cols[r.c1[k]]);
                g.draw(p);
            }
            g.setColor(new Color(170, 170, 190));
            g.setFont(g.getFont().deriveFont(10f));
            g.drawString(String.format(Locale.ROOT, "kernels: distance 0 → %.0f cells", maxD), 4, getHeight() - 2);
        }
    }
}
