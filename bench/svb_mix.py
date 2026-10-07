import mmap, struct, sys
out = sys.argv[1]; parts = [(a.split(":")[0], int(a.split(":")[1]), int(a.split(":")[2])) for a in sys.argv[2:]]
n = total = longest = 0; off = []
with open(out, "wb") as o:
    o.write(b"\0" * 32); pos = 32
    for src, skip, take in parts:
        with open(src, "rb") as f:
            m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            p, seen, got, end = 32, 0, 0, len(m)
            while p < end and got < take:
                tag = struct.unpack_from("<I", m, p)[0]; ln = tag & 0x7fffffff
                sz = 4 + (ln if tag >> 31 else (ln + 3) // 4)
                if seen >= skip:
                    o.write(m[p:p + sz]); off.append(pos); pos += sz; got += ln; total += ln; n += 1; longest = max(longest, ln)
                seen += ln; p += sz
            if got < take: sys.exit(f"{src}: only {got} of {take} bases after skipping {skip}")
    off.append(pos)
    o.seek(0); o.write(b"SVB1" + struct.pack("<IQQ", 1, n, total) + b"\0" * 8)
with open(out + ".idx", "wb") as o:
    o.write(b"SVBI" + struct.pack("<IQQ", 1, n, longest) + b"\0" * 8)
    o.write(struct.pack("<%dQ" % len(off), *off))
print(out, n, "reads", total, "bases")
