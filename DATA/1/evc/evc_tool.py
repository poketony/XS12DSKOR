#!/usr/bin/env python3
"""
evc_tool.py - Xenosaga EVC 파일 텍스트 추출/삽입 툴
사용법:
  python3 evc_tool.py extract <input.evc> [vcf.json]
  python3 evc_tool.py insert  <original.evc> <input.txt> [vcf.json]
"""

import sys, json, struct

# ── Shift-JIS 판별 ───────────────────────────────────────────────────────────

def is_shiftjis_char(b1, b2):
    if 0xA1 <= b1 <= 0xDF: return 1
    if (0x81 <= b1 <= 0x9F) or (0xE0 <= b1 <= 0xFC):
        if (0x40 <= b2 <= 0x7E) or (0x80 <= b2 <= 0xFC): return 2
    if 0x20 <= b1 <= 0x7E: return 1
    return 0

def is_string_jis(data, start, length):
    i = 0
    while i < length:
        b1 = data[start+i] if (start+i) < len(data) else 0
        b2 = data[start+i+1] if (start+i+1) < len(data) else 0
        check = is_shiftjis_char(b1, b2)
        if not check: return False
        i += check
    return True

# ── 파일 파싱 ────────────────────────────────────────────────────────────────

def parse_evc(data):
    entries = []
    jumps = []
    i = 0
    n = len(data)
    while i < n:
        b = data[i]
        # 선택지 블록: 23 [count] 00
        if b == 0x23 and (i+2) < n and data[i+2] == 0x00:
            count = data[i+1]
            j = i+3
            found = 0
            while j < n and data[j] == 0x24 and found < count:
                if j+2 >= n: break
                idx = data[j+1]; slen = data[j+2]
                text_start = j+3; text_end = text_start+slen
                if text_end >= n or data[text_end] != 0x00: break
                entries.append({
                    'type': 'choice', 'offset': j, 'len_offset': j+2,
                    'text_start': j+3, 'idx': idx,
                    'text_bytes': bytes(data[text_start:text_end]),
                })
                j = text_end+1; found += 1
            i = j; continue
        # 점프 opcode: 26 [idx] [4바이트 LE 절대오프셋]
        if b == 0x26 and (i+5) < n:
            idx = data[i+1]
            offset_val = struct.unpack_from('<I', data, i+2)[0]
            if offset_val < n and idx < 16:
                jumps.append({'pos': i, 'idx': idx, 'offset': offset_val})
        if b == 0x27 and (i+4) < n:
            offset_val = struct.unpack_from('<I', data, i+1)[0]
            if offset_val < n:
                jumps.append({'pos': i, 'idx': -1, 'offset': offset_val, 'op': 0x27})
        # 일반 대사: [len] [text] [00]
        if b != 0x00:
            s = b-256 if b > 127 else b
            if s > 0:
                end = i+s+1
                if end < n and data[end] == 0x00 and is_string_jis(data, i+1, s):
                    entries.append({
                        'type': 'dialog', 'offset': i, 'len_offset': i,
                        'text_start': i+1, 'text_bytes': bytes(data[i+1:i+1+s]),
                    })
        i += 1
    return entries, jumps

# ── VCF ─────────────────────────────────────────────────────────────────────

def load_vcf(path):
    with open(path, encoding='utf-8-sig') as f:
        data = json.load(f)
    tbl = data['replace-table']
    kor2sjis, sjis2kor = {}, {}
    for kor, kanji in tbl.items():
        try:
            sb = kanji.encode('shift-jis')
            kor2sjis[kor] = sb; sjis2kor[sb] = kor
        except: pass
    return kor2sjis, sjis2kor

def sjis_to_kor(text_bytes, sjis2kor):
    result = []; i = 0
    while i < len(text_bytes):
        b1 = text_bytes[i]
        if (0x81 <= b1 <= 0x9F) or (0xE0 <= b1 <= 0xFC):
            if i+1 < len(text_bytes):
                key = bytes([b1, text_bytes[i+1]])
                if key in sjis2kor: result.append(sjis2kor[key]); i += 2; continue
                try: result.append(key.decode('shift-jis'))
                except: result.append(f'[{b1:02x}{text_bytes[i+1]:02x}]')
                i += 2; continue
        key = bytes([b1])
        if key in sjis2kor: result.append(sjis2kor[key])
        else:
            try: result.append(bytes([b1]).decode('shift-jis'))
            except: result.append(f'[{b1:02x}]')
        i += 1
    return ''.join(result)

def kor_to_sjis(text_str, kor2sjis):
    result = bytearray()
    for ch in text_str:
        if ch in kor2sjis: result += kor2sjis[ch]
        else:
            try: result += ch.encode('shift-jis')
            except: result += b'?'
    return bytes(result)

# ── 점프 오프셋 보정 ─────────────────────────────────────────────────────────

def fix_jump_offsets(data, insert_at, diff):
    i = 0
    n = len(data)
    while i < n-5:
        if data[i] == 0x26:
            idx = data[i+1]
            offset_val = struct.unpack_from('<I', data, i+2)[0]
            if offset_val < n and idx < 16:
                if offset_val >= insert_at:
                    struct.pack_into('<I', data, i+2, offset_val + diff)
        elif data[i] == 0x27:
            offset_val = struct.unpack_from('<I', data, i+1)[0]
            if 0 < offset_val < n:
                if offset_val >= insert_at:
                    struct.pack_into('<I', data, i+1, offset_val + diff)
        i += 1

# ── EXTRACT ──────────────────────────────────────────────────────────────────

def cmd_extract(evc_path, txt_path, vcf_path=None):
    with open(evc_path, 'rb') as f:
        data = f.read()
    sjis2kor = load_vcf(vcf_path)[1] if vcf_path else {}
    entries, _ = parse_evc(data)
    with open(txt_path, 'w', encoding='utf-8') as out:
        for e in entries:
            raw = e['text_bytes']
            text = sjis_to_kor(raw, sjis2kor) if vcf_path else raw.decode('shift-jis', errors='replace')
            tag = 'C' if e['type'] == 'choice' else 'D'
            idx_str = f":{e['idx']}" if e['type'] == 'choice' else ''
            out.write(f"[{tag}{idx_str}@{e['offset']:08x}]{text}\n")
    print(f"추출 완료: {len(entries)}개 → {txt_path}")

# ── INSERT ───────────────────────────────────────────────────────────────────

def cmd_insert(ori_path, txt_path, out_path, vcf_path=None):
    with open(ori_path, 'rb') as f:
        data = bytearray(f.read())
    kor2sjis = load_vcf(vcf_path)[0] if vcf_path else {}

    replacements = {}
    with open(txt_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.startswith('['): continue
            bracket_end = line.index(']')
            tag_part = line[1:bracket_end]
            text_part = line[bracket_end+1:]
            at_pos = tag_part.index('@')
            offset = int(tag_part[at_pos+1:], 16)
            if vcf_path and kor2sjis:
                new_bytes = kor_to_sjis(text_part, kor2sjis)
            else:
                new_bytes = text_part.encode('shift-jis', errors='replace')
            replacements[offset] = new_bytes

    entries, _ = parse_evc(data)
    to_patch = [e for e in entries if e['offset'] in replacements]
    to_patch.sort(key=lambda e: e['offset'], reverse=True)

    patched = skipped = 0
    for e in to_patch:
        new_text = replacements[e['offset']]
        new_len = len(new_text)
        old_len = len(e['text_bytes'])

        if new_len > 127:
            print(f"  SKIP 0x{e['offset']:08x}: 너무 김 ({new_len}바이트)")
            skipped += 1
            continue

        len_pos    = e['len_offset']
        text_start = e['text_start']
        text_end   = text_start + old_len

        if new_len == old_len:
            data[text_start:text_end] = new_text
        elif new_len < old_len:
            diff = new_len - old_len
            data[len_pos] = new_len
            data = (data[:text_start] + bytearray(new_text) + b'\x00' + data[text_end+1:])
            fix_jump_offsets(data, text_end+1, diff)
        else:
            diff = new_len - old_len
            data[len_pos] = new_len
            data = (data[:text_start] + bytearray(new_text) + b'\x00' + data[text_end+1:])
            fix_jump_offsets(data, text_end+1, diff)

        patched += 1

    with open(out_path, 'wb') as f:
        f.write(data)
    print(f"삽입 완료: {patched}개 패치, {skipped}개 스킵 → {out_path}")
    print(f"파일 크기: {len(data)}바이트")

# ── MAIN ─────────────────────────────────────────────────────────────────────

def find_json_in_dir(base_path):
    import os
    folder = os.path.dirname(os.path.abspath(base_path))
    jsons = [f for f in os.listdir(folder) if f.endswith('.json')]
    if len(jsons) == 1:
        return os.path.join(folder, jsons[0])
    elif len(jsons) > 1:
        print(f"  경고: .json 파일이 여러 개입니다 → {jsons}")
        print(f"  첫 번째 파일 사용: {jsons[0]}")
        return os.path.join(folder, jsons[0])
    return None

def has_korean(text):
    return any('\uAC00' <= ch <= '\uD7A3' for ch in text)

def usage():
    print(__doc__); sys.exit(1)

if __name__ == '__main__':
    import os
    args = sys.argv[1:]
    if len(args) < 2: usage()
    cmd = args[0]

    if cmd == 'extract':
        evc_path = args[1]
        txt_path = os.path.splitext(evc_path)[0] + '.txt'
        vcf = args[2] if len(args) >= 3 and args[2].endswith('.json') else None
        cmd_extract(evc_path, txt_path, vcf)

    elif cmd == 'insert':
        if len(args) < 3: usage()
        evc_path = args[1]
        txt_path = args[2]
        out_path = evc_path + '.new'
        vcf = None
        if len(args) >= 4 and args[3].endswith('.json'):
            vcf = args[3]
        else:
            with open(txt_path, 'r', encoding='utf-8') as f:
                sample = f.read(4096)
            if has_korean(sample):
                vcf = find_json_in_dir(evc_path)
                if vcf:
                    print(f"  한글 감지 → VCF 자동 사용: {os.path.basename(vcf)}")
                else:
                    print("  경고: 한글 감지됐지만 같은 폴더에 .json 파일 없음")
        cmd_insert(evc_path, txt_path, out_path, vcf)

    else:
        usage()