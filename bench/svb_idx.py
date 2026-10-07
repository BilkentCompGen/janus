import mmap, struct, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src, "rb") as f:
    m = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    n = struct.unpack_from("<Q", m, 8)[0]
    pos, longest, off = 32, 0, []
    for _ in range(n):
        tag = struct.unpack_from("<I", m, pos)[0]
        ln = tag & 0x7fffffff
        longest = max(longest, ln)
        off.append(pos)
        pos += 4 + (ln if tag >> 31 else (ln + 3) // 4)
    off.append(pos)
    assert pos == len(m), (pos, len(m))
with open(dst, "wb") as o:
    o.write(b"SVBI" + struct.pack("<IQQ", 1, n, longest) + b"\0" * 8)
    o.write(struct.pack("<%dQ" % len(off), *off))
print(src, n, "reads, longest", longest)
