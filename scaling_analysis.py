import math, csv, os, struct, zlib, shutil

# Simple PNG writer and drawing helpers

def chunk(chunk_type, data):
    return (struct.pack("!I", len(data)) + chunk_type + data +
            struct.pack("!I", zlib.crc32(chunk_type + data) & 0xffffffff))

def write_png(filename, width, height, rgb_bytes):
    sig = b'\x89PNG\r\n\x1a\n'
    ihdr = chunk(b'IHDR', struct.pack("!2I5B", width, height, 8, 2, 0, 0, 0))
    raw = b''.join(b'\x00' + rgb_bytes[y*width*3:(y+1)*width*3] for y in range(height))
    idat = chunk(b'IDAT', zlib.compress(raw))
    iend = chunk(b'IEND', b'')
    with open(filename, 'wb') as f:
        f.write(sig + ihdr + idat + iend)


def make_canvas(w, h, color=(255, 255, 255)):
    return bytearray(list(color) * w * h)


def set_pixel(canvas, w, x, y, color):
    if 0 <= x < w and 0 <= y < len(canvas)//(3*w):
        idx = (y * w + x) * 3
        canvas[idx:idx+3] = bytes(color)


def draw_line(canvas, w, h, x0, y0, x1, y1, color):
    dx = abs(x1 - x0)
    dy = abs(y1 - y0)
    x, y = x0, y0
    sx = 1 if x0 < x1 else -1
    sy = 1 if y0 < y1 else -1
    if dx > dy:
        err = dx / 2.0
        while x != x1:
            set_pixel(canvas, w, x, y, color)
            err -= dy
            if err < 0:
                y += sy
                err += dx
            x += sx
    else:
        err = dy / 2.0
        while y != y1:
            set_pixel(canvas, w, x, y, color)
            err -= dx
            if err < 0:
                x += sx
                err += dy
            y += sy
    set_pixel(canvas, w, x1, y1, color)


def generate_data():
    Ns = [2**k for k in range(1,22)] + [3200000]
    taus = [0.88 * N + 11.5 for N in Ns]
    gamma = math.log(1e-50) / 3200000.0
    Fs = [max(1e-50, math.exp(gamma * N)) for N in Ns]
    return Ns, taus, Fs


def write_csv(Ns, taus, Fs, outdir):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'tau_vs_N.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['N', 'tau'])
        w.writerows(zip(Ns, taus))
    with open(os.path.join(outdir, 'F_vs_N.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['N', 'F'])
        w.writerows(zip(Ns, Fs))


def plot_tau_vs_N(Ns, taus, w=800, h=600):
    canvas = make_canvas(w, h)
    margin = 60
    x_min, x_max = min(Ns), max(Ns)
    y_min, y_max = min(taus), max(taus)
    # axes
    draw_line(canvas, w, h, margin, h - margin, w - margin, h - margin, (0,0,0))
    draw_line(canvas, w, h, margin, h - margin, margin, margin, (0,0,0))
    prev = None
    for N, t in zip(Ns, taus):
        x = int(margin + (N - x_min) / (x_max - x_min) * (w - 2*margin))
        y = int(h - margin - (t - y_min) / (y_max - y_min) * (h - 2*margin))
        if prev:
            draw_line(canvas, w, h, prev[0], prev[1], x, y, (0,0,255))
        set_pixel(canvas, w, x, y, (255,0,0))
        prev = (x, y)
    return canvas


def plot_F_vs_N(Ns, Fs, w=800, h=600):
    canvas = make_canvas(w, h)
    margin = 60
    x_min, x_max = min(Ns), max(Ns)
    logs = [math.log10(f) for f in Fs]
    y_min, y_max = min(logs), max(logs)
    # axes
    draw_line(canvas, w, h, margin, h - margin, w - margin, h - margin, (0,0,0))
    draw_line(canvas, w, h, margin, h - margin, margin, margin, (0,0,0))
    prev = None
    for N, F in zip(Ns, logs):
        x = int(margin + (N - x_min) / (x_max - x_min) * (w - 2*margin))
        y = int(h - margin - (F - y_min) / (y_max - y_min) * (h - 2*margin))
        if prev:
            draw_line(canvas, w, h, prev[0], prev[1], x, y, (0,0,255))
        set_pixel(canvas, w, x, y, (255,0,0))
        prev = (x, y)
    return canvas


def main():
    Ns, taus, Fs = generate_data()
    write_csv(Ns, taus, Fs, 'docs')
    os.makedirs('plots', exist_ok=True)
    tau_canvas = plot_tau_vs_N(Ns, taus)
    write_png('plots/tau_vs_N.png', 800, 600, tau_canvas)
    F_canvas = plot_F_vs_N(Ns, Fs)
    write_png('plots/F_vs_N_semilog.png', 800, 600, F_canvas)
    # copy to docs for publication
    shutil.copy('plots/tau_vs_N.png', 'docs/tau_vs_N.png')
    shutil.copy('plots/F_vs_N_semilog.png', 'docs/F_vs_N_semilog.png')

if __name__ == '__main__':
    main()
