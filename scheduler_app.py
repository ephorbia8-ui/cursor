import tkinter as tk
from tkinter import ttk, messagebox, filedialog, colorchooser, simpledialog
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import hashlib, json, os

# OR-Tools CP-SAT
from ortools.sat.python import cp_model

# PDF export (Arabic-enabled)
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

import arabic_reshaper
from bidi.algorithm import get_display

DAYS = ["Sun","Mon","Tue","Wed","Thu"]
SLOT_START = 8
SLOT_COUNT = 8
SLOTS = list(range(SLOT_START, SLOT_START + SLOT_COUNT))  # 08..15 start-hours
MORNING_CUTOFF = 12
AFTERNOON_START = 12
OUTPUT_DIR = "outputs"
FONTS_DIR = os.path.join(os.path.dirname(__file__), "fonts")
DEFAULT_AR_FONT_NAME = "Amiri"
DEFAULT_AR_FONT_FILE = "Amiri-Regular.ttf"

@dataclass
class CourseInfo:
  name: str
  students: int
  max_group: int
  theory: int
  lab: int
  color: str
  groups: int
  total_hours: int
  hall: Optional[str] = None
  lecturer: Optional[str] = None

@dataclass
class Session:
  course: str
  group: int
  kind: str            # 'theory' | 'lab'
  duration: int        # hours as slots
  day: Optional[str] = None
  start_hour: Optional[int] = None  # one of SLOTS
  priority_ok: bool = True

def deterministic_color(name: str) -> str:
  h = hashlib.sha1(name.encode('utf-8')).hexdigest()
  v = int(h[:6], 16)
  r = (v >> 16) & 0xFF
  g = (v >> 8) & 0xFF
  b = v & 0xFF
  r = (r + 220)//2
  g = (g + 220)//2
  b = (b + 220)//2
  return f"#{r:02x}{g:02x}{b:02x}"

def lighten_color(hexcolor: str, factor: float = 0.6) -> str:
  hexcolor = hexcolor.lstrip("#")
  r = int(hexcolor[0:2], 16)
  g = int(hexcolor[2:4], 16)
  b = int(hexcolor[4:6], 16)
  nr = int(r + (255 - r) * factor)
  ng = int(g + (255 - g) * factor)
  nb = int(b + (255 - b) * factor)
  return f"#{nr:02x}{ng:02x}{nb:02x}"

def ceil_div(a:int,b:int)->int:
  return (a + b - 1)//b if b>0 else 1

def split_theory_blocks(hours:int)->List[int]:
  if hours <= 0:
    return []
  if hours == 2:
    return [2]
  if hours == 4:
    return [2, 2]
  blocks=[]
  rem=hours
  while rem>=2:
    blocks.append(2); rem-=2
  if rem==1:
    blocks.append(1)
  return blocks

def split_lab_block(hours:int)->List[int]:
  return [hours] if hours > 0 else []

def slot_idx_from_hour(h: int) -> int:
  return h - SLOT_START

def slot_label(slot_idx: int) -> str:
  start_m = SLOT_START*60 + slot_idx*60
  end_m = start_m + 45
  sh, sm = start_m//60, start_m%60
  eh, em = end_m//60, end_m%60
  return f"{sh:02d}:{sm:02d}-{eh:02d}:{em:02d}"

def hex_to_rl_color(hex_color: str):
  return colors.HexColor(hex_color)

def ensure_arabic_font(root: tk.Tk) -> str:
  os.makedirs(FONTS_DIR, exist_ok=True)
  font_path = os.path.join(FONTS_DIR, DEFAULT_AR_FONT_FILE)
  if not os.path.exists(font_path):
    try:
      import urllib.request
      url = "https://github.com/aliftype/amiri/releases/download/0.116/Amiri-Regular.ttf"
      urllib.request.urlretrieve(url, font_path)
    except Exception:
      # Do not prompt user; fallback to base fonts (Arabic shaping may degrade)
      return "Helvetica"
  try:
    pdfmetrics.registerFont(TTFont(DEFAULT_AR_FONT_NAME, font_path))
    return DEFAULT_AR_FONT_NAME
  except Exception:
    return "Helvetica"

def ar_text(txt: str) -> str:
  if not txt:
    return ""
  try:
    reshaped = arabic_reshaper.reshape(txt)
    return get_display(reshaped)
  except Exception:
    return txt

class SchedulerEngine:
  def __init__(self, days=DAYS, slots=SLOTS):
    self.days = days[:]
    self.slots = slots[:]
    self.morning_cutoff = MORNING_CUTOFF
    self.afternoon_start = AFTERNOON_START

  def schedule(self, courses: List[CourseInfo], free_day_overrides: Optional[Dict[int, str]] = None):
    max_groups = max((c.groups for c in courses), default=1)
    all_groups = list(range(1, max_groups + 1))

    offday: Dict[int, str] = {}
    for g in all_groups:
      if free_day_overrides and g in free_day_overrides and free_day_overrides[g] in self.days:
        offday[g] = free_day_overrides[g]
      else:
        offday[g] = self.days[(g-1) % len(self.days)]

    # Build blocks per course (scheduled ONCE per course)
    course_blocks: Dict[str, List[Tuple[str,int]]] = {}
    for c in courses:
      blocks: List[Tuple[str,int]] = []
      for b in split_lab_block(c.lab):
        blocks.append(("lab", b))
      for b in split_theory_blocks(c.theory):
        blocks.append(("theory", b))
      course_blocks[c.name] = blocks

    model = cp_model.CpModel()

    # Decision vars: s[c,b,d,h] = 1 if block b of course c starts at (d,h)
    s: Dict[Tuple[str,int,str,int], cp_model.IntVar] = {}
    for c in courses:
      cname = c.name
      blocks = course_blocks[cname]
      for b_idx,(kind,dur) in enumerate(blocks):
        valid = []
        for d in self.days:
          for h in self.slots:
            v = model.NewBoolVar(f"s[{cname},B{b_idx},{d},{h}]")
            s[(cname,b_idx,d,h)] = v
            if h + dur <= self.slots[-1] + 1:
              valid.append(v)
            else:
              model.Add(v == 0)
        # Exactly one start per block
        model.Add(sum(valid) == 1)

    # No-overlap across courses: at any (d,h_unit) at most one running block
    for d in self.days:
      for h_unit in self.slots:
        running = []
        for c in courses:
          cname = c.name
          for b_idx,(kind,dur) in enumerate(course_blocks[cname]):
            for h in self.slots:
              if h <= h_unit <= h + dur - 1 and h + dur <= self.slots[-1] + 1:
                running.append(s[(cname,b_idx,d,h)])
        if running:
          model.Add(sum(running) <= 1)

    # Preferences: theory in morning, lab in afternoon; off-day as soft penalty per group
    penalties = []
    for c in courses:
      cname = c.name
      for b_idx,(kind,dur) in enumerate(course_blocks[cname]):
        for d in self.days:
          for h in self.slots:
            if h + dur <= self.slots[-1] + 1:
              v = model.NewIntVar(0, 1, f"pref[{cname},B{b_idx},{d},{h}]")
              model.Add(v == 1).OnlyEnforceIf(s[(cname,b_idx,d,h)])
              model.Add(v == 0).OnlyEnforceIf(s[(cname,b_idx,d,h)].Not())
              if kind == "theory":
                if h >= MORNING_CUTOFF:
                  penalties.append(v)
              else:
                if h < AFTERNOON_START:
                  penalties.append(v)
              # Off-day penalty for each group that has this day off
              for g in all_groups:
                if d == offday[g]:
                  penalties.append(v)

    if penalties:
      model.Minimize(sum(penalties))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = 20.0
    solver.parameters.num_search_workers = 8
    status = solver.Solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
      return [], {}

    # Build identical schedule for all groups
    per_group: Dict[int, List[Session]] = {g: [] for g in all_groups}
    for c in courses:
      cname = c.name
      for b_idx,(kind,dur) in enumerate(course_blocks[cname]):
        for d in self.days:
          for h in self.slots:
            if h + dur <= self.slots[-1] + 1 and solver.BooleanValue(s[(cname,b_idx,d,h)]):
              for g in all_groups:
                per_group[g].append(Session(course=cname, group=g, kind=kind, duration=dur, day=d, start_hour=h, priority_ok=True))

    def sort_key(ses: Session):
      return (DAYS.index(ses.day) if ses.day in DAYS else 999, ses.start_hour or 999, ses.course, ses.kind)
    for g in per_group:
      per_group[g].sort(key=sort_key)
    all_sessions = []
    for g in per_group:
      all_sessions.extend(per_group[g])
    all_sessions.sort(key=sort_key)
    return all_sessions, per_group

class SchedulerGUI:
  def __init__(self, root):
    self.root = root
    root.title("Advanced Scheduler (Sun-Thu, 8 periods)")
    root.geometry("1220x780")

    self.courses: Dict[str, CourseInfo] = {}
    self.engine = SchedulerEngine(days=DAYS, slots=SLOTS)
    self.last_scheduled: List[Session] = []
    self.last_per_group: Dict[int, List[Session]] = {}
    self.multislot_groups: Dict[Tuple[str,str,int,int], List[int]] = {}

    self.selected_color = None
    self.free_day_overrides: Dict[int, str] = {}

    self.temp_hall: Optional[str] = None
    self.temp_lecturer: Optional[str] = None

    self.hall_label_var = tk.StringVar(value="—")
    self.lect_label_var = tk.StringVar(value="—")

    self.build_ui()

  def build_ui(self):
    top = ttk.Frame(self.root, padding=8)
    top.pack(side="top", fill="x")

    ttk.Label(top, text="Course name").grid(row=0,column=0, sticky="w")
    self.e_name = ttk.Entry(top, width=20); self.e_name.grid(row=0,column=1, padx=4)

    ttk.Label(top, text="Students").grid(row=0,column=2, sticky="w")
    self.e_students = ttk.Entry(top, width=8); self.e_students.grid(row=0,column=3, padx=4)

    ttk.Label(top, text="Max per group").grid(row=0,column=4, sticky="w")
    self.e_max = ttk.Entry(top, width=8); self.e_max.grid(row=0,column=5, padx=4)

    ttk.Label(top, text="Theory hrs").grid(row=1,column=0, sticky="w")
    self.e_th = ttk.Entry(top, width=8); self.e_th.grid(row=1,column=1, padx=4)

    ttk.Label(top, text="Lab hrs").grid(row=1,column=2, sticky="w")
    self.e_lab = ttk.Entry(top, width=8); self.e_lab.grid(row=1,column=3, padx=4)

    ttk.Button(top, text="Pick color", command=self.pick_color).grid(row=0,column=6, padx=4)
    self.color_canvas = tk.Canvas(top, width=36, height=18, bg="#ffffff", bd=1, relief="sunken"); self.color_canvas.grid(row=0,column=7)

    ttk.Button(top, text="قاعة المحاضرات", command=self.set_hall).grid(row=1,column=6, padx=4, sticky="w")
    ttk.Label(top, textvariable=self.hall_label_var).grid(row=1,column=7, padx=4, sticky="w")

    ttk.Button(top, text="المحاضر", command=self.set_lecturer).grid(row=1,column=8, padx=4, sticky="w")
    ttk.Label(top, textvariable=self.lect_label_var).grid(row=1,column=9, padx=4, sticky="w")

    ttk.Button(top, text="Add course", command=self.add_course).grid(row=0,column=8, padx=6)
    ttk.Button(top, text="Edit selected", command=self.edit_selected).grid(row=0,column=9, padx=6)
    ttk.Button(top, text="Delete selected", command=self.delete_selected).grid(row=0,column=10, padx=6)

    ttk.Button(top, text="Save list", command=self.save_courses).grid(row=2,column=0, pady=6)
    ttk.Button(top, text="Load list", command=self.load_courses).grid(row=2,column=1, pady=6)
    ttk.Button(top, text="Free day settings", command=self.open_free_day_settings).grid(row=2,column=3, pady=6)
    ttk.Button(top, text="Generate schedule", command=self.generate).grid(row=2,column=4, pady=6)
    ttk.Button(top, text="Export PDFs", command=self.export_dialog).grid(row=2,column=5, pady=6)

    # Improve spacing in notebook cells: larger default sizes
    self.cell_width = 170
    self.cell_height = 56

    mid = ttk.Frame(self.root, padding=8); mid.pack(fill="x")
    cols = ("students","maxg","groups","th","lab","total")
    self.tree = ttk.Treeview(mid, columns=cols, show="headings", height=8)
    for k,h in zip(cols, ("Students","MaxGrp","#Groups","Th","Lab","TotalHrs")):
      self.tree.heading(k, text=h)
    self.tree.pack(side="left", fill="x", expand=True)
    self.tree.bind("<<TreeviewSelect>>", self.on_select)
    sb = ttk.Scrollbar(mid, orient="vertical", command=self.tree.yview); sb.pack(side="right", fill="y"); self.tree.configure(yscroll=sb.set)

    self.notebook = ttk.Notebook(self.root); self.notebook.pack(fill="both", expand=True, padx=8, pady=8)
    self.status = ttk.Label(self.root, text="Ready"); self.status.pack(side="bottom", fill="x")

  def pick_color(self):
    v = colorchooser.askcolor()
    if v and v[1]:
      self.selected_color = v[1]
      self.color_canvas.config(bg=self.selected_color)

  def set_hall(self):
    val = simpledialog.askstring("قاعة المحاضرات", "ادخل اسم القاعة:", parent=self.root)
    if val is not None:
      self.temp_hall = val.strip() or None
      self.hall_label_var.set(self.temp_hall or "—")

  def set_lecturer(self):
    val = simpledialog.askstring("المحاضر", "ادخل اسم المحاضر:", parent=self.root)
    if val is not None:
      self.temp_lecturer = val.strip() or None
      self.lect_label_var.set(self.temp_lecturer or "—")

  def add_course(self):
    name = self.e_name.get().strip()
    try:
      students = int(self.e_students.get().strip())
      maxg = int(self.e_max.get().strip())
      th = int(self.e_th.get().strip() or 0)
      lab = int(self.e_lab.get().strip() or 0)
    except Exception:
      messagebox.showerror("Input error","Enter integers for students/max/hour fields")
      return
    if not name:
      messagebox.showerror("Input error","Course name required"); return
    if maxg <= 0:
      messagebox.showerror("Input error","Max per group must be > 0"); return
    if any(x < 0 for x in (students, th, lab)):
      messagebox.showerror("Input error","Negative values not allowed"); return

    color = self.selected_color or deterministic_color(name)
    groups = ceil_div(students, max(1, maxg))
    total = groups * (th + lab)
    self.courses[name] = CourseInfo(name=name, students=students, max_group=maxg, theory=th, lab=lab, color=color, groups=groups, total_hours=total, hall=self.temp_hall, lecturer=self.temp_lecturer)
    self.refresh_tree()
    self.clear_inputs()

  def clear_inputs(self):
    self.e_name.delete(0, tk.END); self.e_students.delete(0, tk.END); self.e_max.delete(0, tk.END)
    self.e_th.delete(0, tk.END); self.e_lab.delete(0, tk.END)
    self.selected_color = None; self.color_canvas.config(bg="#ffffff")
    self.temp_hall = None; self.hall_label_var.set("—")
    self.temp_lecturer = None; self.lect_label_var.set("—")

  def refresh_tree(self):
    for r in self.tree.get_children():
      self.tree.delete(r)
    for nm,info in self.courses.items():
      self.tree.insert("", tk.END, iid=nm, values=(info.students, info.max_group, info.groups, info.theory, info.lab, info.total_hours))

  def on_select(self, event=None):
    sel = self.tree.selection()
    if not sel: return
    iid = sel[0]
    info = self.courses.get(iid)
    if info:
      self.e_name.delete(0, tk.END); self.e_name.insert(0, info.name)
      self.e_students.delete(0, tk.END); self.e_students.insert(0, str(info.students))
      self.e_max.delete(0, tk.END); self.e_max.insert(0, str(info.max_group))
      self.e_th.delete(0, tk.END); self.e_th.insert(0, str(info.theory))
      self.e_lab.delete(0, tk.END); self.e_lab.insert(0, str(info.lab))
      self.selected_color = info.color; self.color_canvas.config(bg=info.color)
      self.temp_hall = info.hall; self.hall_label_var.set(info.hall or "—")
      self.temp_lecturer = info.lecturer; self.lect_label_var.set(info.lecturer or "—")

  def edit_selected(self):
    sel = self.tree.selection()
    if not sel:
      messagebox.showinfo("Edit","Select a course"); return
    old_name = sel[0]
    name = self.e_name.get().strip()
    try:
      students = int(self.e_students.get().strip()); maxg = int(self.e_max.get().strip())
      th = int(self.e_th.get().strip() or 0); lab = int(self.e_lab.get().strip() or 0)
    except Exception:
      messagebox.showerror("Input error","Integers only"); return
    if not name:
      messagebox.showerror("Input error","Course name required"); return
    if maxg <= 0:
      messagebox.showerror("Input error","Max per group must be > 0"); return
    if name != old_name and name in self.courses:
      messagebox.showerror("Edit", "Course with this name exists"); return
    color = self.selected_color or deterministic_color(name)
    groups = ceil_div(students, max(1, maxg)); total = groups * (th + lab)
    del self.courses[old_name]
    self.courses[name] = CourseInfo(name=name, students=students, max_group=maxg, theory=th, lab=lab, color=color, groups=groups, total_hours=total, hall=self.temp_hall, lecturer=self.temp_lecturer)
    self.refresh_tree(); self.clear_inputs()

  def delete_selected(self):
    sel = self.tree.selection()
    if not sel:
      messagebox.showinfo("Delete","Select a course"); return
    iid = sel[0]; del self.courses[iid]; self.refresh_tree()

  def save_courses(self):
    if not self.courses:
      messagebox.showinfo("Save","No courses to save"); return
    filename = filedialog.asksaveasfilename(defaultextension=".json", filetypes=[("JSON files", "*.json")])
    if filename:
      data = []
      for info in self.courses.values():
        data.append({
          "name": info.name, "students": info.students, "max_group": info.max_group,
          "theory": info.theory, "lab": info.lab, "color": info.color,
          "groups": info.groups, "total_hours": info.total_hours,
          "hall": info.hall, "lecturer": info.lecturer
        })
      with open(filename, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
      messagebox.showinfo("Save", f"Saved {len(data)} courses to {filename}")

  def load_courses(self):
    filename = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
    if filename:
      try:
        with open(filename, "r", encoding="utf-8") as f:
          data = json.load(f)
        self.courses.clear()
        for item in data:
          info = CourseInfo(
            name=item["name"], students=item["students"], max_group=item["max_group"],
            theory=item["theory"], lab=item["lab"], color=item["color"],
            groups=item["groups"], total_hours=item["total_hours"],
            hall=item.get("hall"), lecturer=item.get("lecturer")
          )
          self.courses[info.name] = info
        self.refresh_tree()
        messagebox.showinfo("Load", f"Loaded {len(data)} courses from {filename}")
      except Exception as e:
        messagebox.showerror("Load error", f"Failed to load: {e}")

  def open_free_day_settings(self):
    if not self.courses:
      messagebox.showinfo("Free day","Add courses first"); return
    max_groups = max((c.groups for c in self.courses.values()), default=1)
    win = tk.Toplevel(self.root); win.title("Free day settings")
    rows = ttk.Frame(win, padding=8); rows.pack(fill="both", expand=True)
    combos: Dict[int, ttk.Combobox] = {}
    for g in range(1, max_groups+1):
      ttk.Label(rows, text=f"Group {g}").grid(row=g, column=0, padx=6, pady=4, sticky="w")
      cb = ttk.Combobox(rows, values=DAYS, state="readonly", width=8)
      cb.set(DAYS[(g-1)%len(DAYS)])
      cb.grid(row=g, column=1, padx=6, pady=4, sticky="w")
      combos[g] = cb
    def apply():
      self.free_day_overrides = {g: combos[g].get() for g in combos}
      win.destroy()
    ttk.Button(rows, text="Apply", command=apply).grid(row=max_groups+1, column=0, columnspan=2, pady=8)

  def generate(self):
    if not self.courses:
      messagebox.showinfo("Generate","Add courses first"); return
    all_sessions, per_group = self.engine.schedule(list(self.courses.values()), free_day_overrides=self.free_day_overrides)
    if not all_sessions:
      self.status.config(text="Solver found no feasible solution. Try adjusting inputs.")
      return
    self.last_scheduled = all_sessions
    self.last_per_group = per_group
    multi: Dict[Tuple[str,str,int,int], List[int]] = {}
    for s in all_sessions:
      key = (s.course, s.day, s.start_hour, s.duration)
      multi.setdefault(key, []).append(s.group)
    for k in multi:
      multi[k] = sorted(set(multi[k]))
    self.multislot_groups = multi
    self.render_calendar(per_group)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    try:
      self.write_pdfs(self.last_scheduled, self.last_per_group, OUTPUT_DIR)
      self.status.config(text=f"Generated. PDFs in '{OUTPUT_DIR}/'")
    except Exception as e:
      self.status.config(text=f"PDF export failed: {e}")

  def render_calendar(self, per_group: Dict[int, List[Session]]):
    for t in self.notebook.tabs():
      self.notebook.forget(t)
    if not per_group:
      return
    maxg = max(per_group.keys())
    for g in range(1, maxg+1):
      tab = ttk.Frame(self.notebook); self.notebook.add(tab, text=f"Group {g}")
      tk.Label(tab, text="Time\\Day", borderwidth=1, relief="solid", width=12).grid(row=0, column=0, sticky="nsew")
      for j,d in enumerate(DAYS, start=1):
        tk.Label(tab, text=d, borderwidth=1, relief="solid", width=16).grid(row=0, column=j, sticky="nsew")
      cell_map = {}
      for i,h in enumerate(SLOTS, start=1):
        sidx = i-1
        tk.Label(tab, text=slot_label(sidx), borderwidth=1, relief="solid", width=16).grid(row=i, column=0, sticky="nsew")
        for j,d in enumerate(DAYS, start=1):
          frm = tk.Frame(tab, borderwidth=1, relief="solid", width=self.cell_width, height=self.cell_height)
          frm.grid_propagate(False); frm.grid(row=i, column=j, sticky="nsew", padx=1, pady=1)
          cell_map[(d,h)] = frm

      items = per_group.get(g, [])
      drawn = set()
      for s in items:
        if s.day is None or s.start_hour is None:
          continue
        key = (s.course, s.day, s.start_hour, s.duration)
        if key in drawn:
          continue
        drawn.add(key)
        att = self.multislot_groups.get(key, [s.group])
        att_txt = ",".join(f"G{x}" for x in att)
        col = deterministic_color(s.course)
        info = self.courses.get(s.course)
        hall = info.hall if info else None
        lecturer = info.lecturer if info else None
        for k in range(s.duration):
          hh = s.start_hour + k
          frame = cell_map.get((s.day, hh))
          if not frame: continue
          for w in frame.winfo_children(): w.destroy()
          bg = col if s.priority_ok else lighten_color(col, 0.6)
          text = f"{s.course}\n{att_txt}\n{s.kind}"
          extra = []
          if hall: extra.append(hall)
          if lecturer: extra.append(lecturer)
          if extra:
            text += "\n" + " | ".join(extra)
          lbl = tk.Label(frame, text=(text if k==0 else "(cont)"), bg=bg, wraplength=self.cell_width-10, justify="center")
          lbl.pack(fill="both", expand=True)

      for col in range(len(DAYS)+1):
        tab.grid_columnconfigure(col, weight=1)
      for row in range(len(SLOTS)+1):
        tab.grid_rowconfigure(row, weight=1)

  def _build_calendar_table(self, story: List, title: str, items: List[Session], multi: Dict[Tuple[str,str,int,int], List[int]], font_name: str):
    story.append(Paragraph(ar_text(title), ParagraphStyle(name='Title', fontName=font_name, fontSize=16, leading=20, alignment=1)))
    story.append(Spacer(1, 6))

    n_cols = 1 + len(DAYS)
    n_rows = 1 + len(SLOTS)
    data = [["" for _ in range(n_cols)] for __ in range(n_rows)]

    data[0][0] = ar_text("الوقت/اليوم")
    for j,d in enumerate(DAYS, start=1):
      data[0][j] = d
    for i in range(1, n_rows):
      data[i][0] = slot_label(i-1)

    style = TableStyle([
      ('FONT', (0,0), (-1,-1), font_name),
      ('ALIGN', (0,0), (-1,0), 'CENTER'),
      ('ALIGN', (0,1), (0,-1), 'CENTER'),
      ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
      ('GRID', (0,0), (-1,-1), 0.5, colors.black),
      ('BACKGROUND', (0,0), (-1,0), colors.whitesmoke),
    ])

    spans: List[Tuple[Tuple[int,int], Tuple[int,int]]] = []
    for s in items:
      if s.day is None or s.start_hour is None:
        continue
      col = 1 + DAYS.index(s.day)
      row = 1 + slot_idx_from_hour(s.start_hour)
      text = f"{s.course}\n{','.join('G'+str(g) for g in multi.get((s.course, s.day, s.start_hour, s.duration), [s.group]))}\n{s.kind}"
      info = self.courses.get(s.course)
      text_extra = []
      if info and info.hall:
        text_extra.append(info.hall)
      if info and info.lecturer:
        text_extra.append(info.lecturer)
      if text_extra:
        text += "\n" + " | ".join(text_extra)
      data[row][col] = ar_text(text)
      if s.duration > 1:
        spans.append(((col, row), (col, row + s.duration - 1)))
      base_col = deterministic_color(s.course)
      if not s.priority_ok:
        base_col = lighten_color(base_col, 0.6)
      style.add('BACKGROUND', (col, row), (col, row + s.duration - 1), hex_to_rl_color(base_col))

    for (c0,r0),(c1,r1) in spans:
      style.add('SPAN', (c0,r0), (c1,r1))

    total_width = 780
    time_w = 110
    day_w = (total_width - time_w) / len(DAYS)
    col_widths = [time_w] + [day_w]*len(DAYS)

    tbl = Table(data, colWidths=col_widths)
    tbl.setStyle(style)
    story.append(tbl)

  def write_pdfs(self, all_sessions: List[Session], per_group: Dict[int, List[Session]], outdir: str):
    font_name = ensure_arabic_font(self.root)

    multi: Dict[Tuple[str,str,int,int], List[int]] = {}
    for s in all_sessions:
      key = (s.course, s.day, s.start_hour, s.duration)
      multi.setdefault(key, []).append(s.group)
    for k in list(multi.keys()):
      multi[k] = sorted(set(multi[k]))

    master_path = os.path.join(outdir, "master.pdf")
    doc = SimpleDocTemplate(master_path, pagesize=landscape(A4), leftMargin=18, rightMargin=18, topMargin=18, bottomMargin=18)
    story: List = []
    self._build_calendar_table(story, "الجدول الرئيسي", all_sessions, multi, font_name)
    doc.build(story)

    for g, items in per_group.items():
      path = os.path.join(outdir, f"group_{g}.pdf")
      docg = SimpleDocTemplate(path, pagesize=landscape(A4), leftMargin=18, rightMargin=18, topMargin=18, bottomMargin=18)
      story_g: List = []
      self._build_calendar_table(story_g, f"جدول المجموعة {g}", items, multi, font_name)
      total_hours = sum(s.duration for s in items)
      story_g.append(Spacer(1, 8))
      story_g.append(Paragraph(ar_text(f"إجمالي الساعات: {total_hours}"), ParagraphStyle(name='Body', fontName=font_name, fontSize=12, leading=16, alignment=2)))
      docg.build(story_g)

    # All groups in one PDF (each group on its own page)
    all_groups_path = os.path.join(outdir, "all_groups.pdf")
    doc_all = SimpleDocTemplate(all_groups_path, pagesize=landscape(A4), leftMargin=18, rightMargin=18, topMargin=18, bottomMargin=18)
    story_all: List = []
    first = True
    for g, items in per_group.items():
      if not first:
        story_all.append(Spacer(1, 18))
      first = False
      self._build_calendar_table(story_all, f"جدول المجموعة {g}", items, multi, font_name)
      total_hours = sum(s.duration for s in items)
      story_all.append(Spacer(1, 8))
      story_all.append(Paragraph(ar_text(f"إجمالي الساعات: {total_hours}"), ParagraphStyle(name='Body', fontName=font_name, fontSize=12, leading=16, alignment=2)))
    doc_all.build(story_all)

  def export_dialog(self):
    if not self.last_scheduled:
      messagebox.showinfo("Export","Generate schedule first"); return
    messagebox.showinfo("Export", f"PDF files saved in '{OUTPUT_DIR}/' directory")


def main():
  root = tk.Tk()
  app = SchedulerGUI(root)
  root.mainloop()

if __name__ == "__main__":
  main()