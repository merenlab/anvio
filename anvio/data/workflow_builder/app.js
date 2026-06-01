const DATABASE = window.ANVIO_WORKFLOW_DB;
const PARAMETER_RELATIONS = window.ANVIO_PARAMETER_RELATIONS || { versions: {} };
const LOCAL_RUNTIME = window.ANVIO_LOCAL_RUNTIME || null;
const WORLD_WIDTH = 5600;
const WORLD_HEIGHT = 5200;
const NODE_DEFAULT_WIDTH = 230;
const NODE_DEFAULT_HEIGHT = 132;
const DATA_NODE_DEFAULT_HEIGHT = 112;
const HISTORY_LIMIT = 80;
const MIN_ZOOM = 0.45;
const MAX_ZOOM = 1.8;
const ZOOM_IN_STEP = 0.06;
const ZOOM_OUT_STEP = 0.12;
const GITHUB_URL = "https://github.com/merenlab/anvio";

const dom = {
  sourceStatus: document.getElementById("sourceStatus"),
  searchInput: document.getElementById("searchInput"),
  resultList: document.getElementById("resultList"),
  libraryMeta: document.getElementById("libraryMeta"),
  workspace: document.getElementById("workspace"),
  libraryPanel: document.getElementById("libraryPanel"),
  tabs: document.getElementById("tabs"),
  zoomLabel: document.getElementById("zoomLabel"),
  canvasFrame: document.getElementById("canvasFrame"),
  canvasWorld: document.getElementById("canvasWorld"),
  edgeLayer: document.getElementById("edgeLayer"),
  nodeLayer: document.getElementById("nodeLayer"),
  intentBuilder: document.getElementById("intentBuilder"),
  intentBuilderInput: document.getElementById("intentBuilderInput"),
  intentBuilderPreview: document.getElementById("intentBuilderPreview"),
  inspector: document.getElementById("inspector"),
  inspectorPanel: document.getElementById("inspectorPanel"),
  modal: document.getElementById("modal"),
  modalContent: document.getElementById("modalContent"),
  workspaceFileInput: document.getElementById("workspaceFileInput"),
  workflowFileInput: document.getElementById("workflowFileInput")
};

const state = {
  db: DATABASE,
  versionId: "main",
  searchKind: "all",
  searchQuery: "",
  tabs: [],
  activeTabId: "",
  selectedNodeId: "",
  selectedEdgeId: "",
  focusedArtifact: null,
  hoveredEdgeId: "",
  connectStartId: "",
  mode: "select",
  mousePos: { x: 0, y: 0 },
  history: {
    past: [],
    future: []
  },
  settings: {
    fontScale: 1,
    nightMode: false,
    zoom: 1,
    libraryCollapsed: false,
    inspectorCollapsed: false
  }
};

let idCounter = 1;
let textEditSnapshotOpen = false;
let intentPreviewTimer = 0;
const intentGraphCache = new Map();
let importWorkflowDraft = {
  fileName: "",
  text: "",
  selectedFormat: "auto",
  detectedFormat: "unknown",
  commands: [],
  records: [],
  variables: {}
};

function uid(prefix) {
  idCounter += 1;
  return `${prefix}-${Date.now().toString(36)}-${idCounter.toString(36)}`;
}

function currentVersion() {
  return state.db.versions.find((version) => version.id === state.versionId) || state.db.versions[0];
}

function localRuntimeAvailable() {
  return Boolean(LOCAL_RUNTIME?.programs?.length || LOCAL_RUNTIME?.artifacts?.length || LOCAL_RUNTIME?.parameterRelations);
}

function applyLocalRuntimeDatabase(database, runtime) {
  if (!database?.versions?.length || !runtime) return database;
  const baseVersion = localRuntimeBaseVersion(database, runtime);
  const localDocVersion = localRuntimeDocVersion(runtime, baseVersion);
  const localVersion = {
    id: "local",
    label: runtime.version?.label || `local anvi'o${runtime.version?.version ? ` ${runtime.version.version}` : ""}`,
    docsUrl: officialHelpBaseUrl(localDocVersion),
    notice: `Loaded from local anvi'o install${runtime.generatedAt ? ` at ${runtime.generatedAt}` : ""}.`,
    programs: withCompatibleArtifactAliases(mergeLocalPrograms(baseVersion?.programs || [], runtime.programs || [], localDocVersion)),
    dataBlocks: mergeLocalArtifacts(baseVersion?.dataBlocks || [], runtime.artifacts || [], runtime.programs || [], localDocVersion),
    workflows: mergeLocalWorkflows(baseVersion?.workflows || [], runtime.workflows || [], localDocVersion),
    programNetwork: runtime.programNetwork || null,
    localClasses: runtime.classes || [],
    upstreamVersionId: baseVersion?.id || localDocVersion
  };
  const upstreamRelations = PARAMETER_RELATIONS.versions?.[baseVersion?.id] || { programs: {} };
  const localRelations = runtime.parameterRelations || { programs: {} };
  PARAMETER_RELATIONS.versions = PARAMETER_RELATIONS.versions || {};
  PARAMETER_RELATIONS.versions.local = {
    ...upstreamRelations,
    ...localRelations,
    label: localVersion.label,
    programs: {
      ...(upstreamRelations.programs || {}),
      ...(localRelations.programs || {})
    }
  };
  return {
    ...database,
    versions: [localVersion],
    localRuntime: runtime
  };
}

function localRuntimeBaseVersion(database, runtime) {
  const version = String(runtime?.version?.version || "").trim();
  if (version) {
    const exact = database.versions.find((candidate) => String(candidate.id) === version);
    if (exact) return exact;
    const labelMatch = database.versions.find((candidate) => String(candidate.label || "").startsWith(`${version} `));
    if (labelMatch) return labelMatch;
  }
  return database.versions.find((candidate) => candidate.id === "main") || database.versions[0];
}

function localRuntimeDocVersion(runtime, baseVersion) {
  const version = String(runtime?.version?.version || baseVersion?.id || "main").trim();
  return version || "main";
}

function officialHelpBaseUrl(versionId = "main") {
  return `https://anvio.org/help/${encodeURIComponent(versionId || "main")}/`;
}

function officialDocsUrlForItem(kind, id, versionId = "main") {
  const cleanId = String(id || "").trim();
  if (!cleanId) return officialHelpBaseUrl(versionId);
  const section = kind === "data" || kind === "artifact" ? "artifacts" : kind === "workflow" ? "workflows" : "programs";
  return `${officialHelpBaseUrl(versionId)}${section}/${encodeURIComponent(cleanId)}/`;
}

function mergeLocalPrograms(basePrograms, localPrograms, docVersion = "main") {
  const byCommand = new Map(basePrograms.map((program) => [program.command || program.id, clone(program)]));
  localPrograms.forEach((localProgram) => {
    const command = localProgram.command || localProgram.id;
    if (!command) return;
    const existing = byCommand.get(command);
    if (existing) {
      existing.parameters = mergeParameters(existing.parameters || [], localProgram.parameters || []);
      existing.description = localProgram.description || existing.description;
      existing.documentation = localProgram.documentation || existing.documentation;
      existing.docsUrl = localProgram.docsUrl || existing.docsUrl || officialDocsUrlForItem("program", command, docVersion);
      existing.localModule = localProgram.module || localProgram.modulePath || "";
      existing.localSource = true;
      byCommand.set(command, existing);
      return;
    }
    byCommand.set(command, {
      id: command,
      kind: "program",
      title: command,
      command,
      description: localProgram.description || "Program discovered from the local anvi'o installation.",
      requires: localProgram.requires || [],
      provides: localProgram.provides || [],
      parameters: localProgram.parameters || [],
      examples: [command],
      docsUrl: localProgram.docsUrl || officialDocsUrlForItem("program", command, docVersion),
      documentation: localProgram.documentation || "",
      localModule: localProgram.module || localProgram.modulePath || "",
      localSource: true
    });
  });
  return [...byCommand.values()].sort((a, b) => (a.command || a.id).localeCompare(b.command || b.id));
}

function mergeParameters(baseParameters, localParameters) {
  const byName = new Map((baseParameters || []).map((param) => [normalizeParameterName(param.name), clone(param)]));
  localParameters.forEach((localParam) => {
    const key = normalizeParameterName(localParam.name);
    if (!key) return;
    const existing = byName.get(key);
    byName.set(key, existing ? {
      ...existing,
      ...localParam,
      documentation: localParam.documentation || localParam.help || existing.documentation || existing.help || "",
      enabled: Boolean(existing.enabled || localParam.required),
      required: Boolean(existing.required || localParam.required),
      value: existing.value || ""
    } : {
      ...localParam,
      documentation: localParam.documentation || localParam.help || "",
      enabled: Boolean(localParam.required),
      value: ""
    });
  });
  return [...byName.values()].sort((a, b) => normalizeParameterName(a.name).localeCompare(normalizeParameterName(b.name)));
}

function mergeLocalArtifacts(baseArtifacts, localArtifacts, localPrograms = [], docVersion = "main") {
  const byId = new Map(baseArtifacts.map((artifact) => [artifact.id, clone(artifact)]));
  const programIds = new Set([
    ...localPrograms.map((program) => program.command || program.id),
    ...baseArtifacts.filter((artifact) => artifact.kind === "program").map((artifact) => artifact.command || artifact.id)
  ].filter(Boolean));
  localArtifacts.forEach((artifact) => {
    if (!artifact?.id) return;
    if (programIds.has(artifact.id) && /inferred from the local anvi'o installation/i.test(artifact.description || "")) return;
    byId.set(artifact.id, {
      ...(byId.get(artifact.id) || {
        kind: "data",
        title: artifact.id,
        requires: [],
        provides: [],
        parameters: []
      }),
      ...artifact,
      kind: "data",
      docsUrl: artifact.docsUrl || byId.get(artifact.id)?.docsUrl || officialDocsUrlForItem("artifact", artifact.id, docVersion),
      localSource: true
    });
  });
  return [...byId.values()].sort((a, b) => a.id.localeCompare(b.id));
}

function mergeLocalWorkflows(baseWorkflows, localWorkflows, docVersion = "main") {
  const byId = new Map((baseWorkflows || []).map((workflow) => [workflow.id, clone(workflow)]));
  (localWorkflows || []).forEach((workflow) => {
    if (!workflow?.id) return;
    byId.set(workflow.id, {
      ...(byId.get(workflow.id) || {}),
      ...workflow,
      docsUrl: workflow.docsUrl || byId.get(workflow.id)?.docsUrl || officialDocsUrlForItem("workflow", workflow.id, docVersion),
      localSource: true
    });
  });
  return [...byId.values()].sort((a, b) => (a.title || a.id).localeCompare(b.title || b.id));
}

function withCompatibleArtifactAliases(programs) {
  return programs.map((program) => {
    const provides = [...new Set(program.provides || [])];
    if (provides.includes("single-profile-db") && !provides.includes("profile-db")) {
      provides.push("profile-db");
    }
    return { ...program, provides };
  });
}

function getVersionLabel(versionId) {
  return state.db.versions.find((version) => version.id === versionId)?.label || versionId;
}

function currentTab() {
  return state.tabs.find((tab) => tab.id === state.activeTabId) || state.tabs[0];
}

function libraryItems() {
  const version = currentVersion();
  return [
    ...version.programs,
    ...version.workflows,
    ...version.dataBlocks
  ];
}

function findLibraryItem(kind, refId) {
  const version = currentVersion();
  if (kind === "program") return version.programs.find((item) => item.id === refId);
  if (kind === "workflow") return version.workflows.find((item) => item.id === refId);
  return version.dataBlocks.find((item) => item.id === refId) || findDataBlock(refId);
}

function findDataBlock(refId) {
  return currentVersion().dataBlocks.find((item) => item.id === refId) || {
    id: refId,
    kind: "data",
    title: refId,
    artifactType: "DATA",
    description: `${refId} data block`,
    requires: [],
    provides: [],
    parameters: [],
    docsUrl: officialDocsUrlForItem("artifact", refId, currentVersion().upstreamVersionId || currentVersion().id || "main")
  };
}

function activeNodes() {
  return currentTab()?.nodes || [];
}

function activeEdges() {
  return currentTab()?.edges || [];
}

function clone(value) {
  return JSON.parse(JSON.stringify(value));
}

function escapeHtml(value) {
  return String(value ?? "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#39;");
}

function compactText(value, length = 170) {
  const text = String(value || "").replace(/\s+/g, " ").trim();
  if (text.length <= length) return text;
  return `${text.slice(0, length - 1).trim()}...`;
}

function safeName(value) {
  return String(value || "item")
    .replace(/^--?/, "")
    .replace(/[^A-Za-z0-9_]+/g, "_")
    .replace(/^_+|_+$/g, "")
    .toLowerCase() || "item";
}

function cleanDocText(text) {
  if (!text) return "";
  const noise = "Mentioned in the official documentation; the official help page does not publish a complete argparse/default table.";
  return text === noise ? "" : text;
}

function normalizeParameterName(value) {
  const text = String(value || "").trim().toLowerCase();
  if (!text) return "";
  if (text.startsWith("-")) return text;
  return `--${text.replace(/^_+/, "").replaceAll("_", "-")}`;
}

function parameterNameCandidates(param) {
  return [
    param?.name,
    ...(param?.aliases || [])
  ].map(normalizeParameterName).filter(Boolean);
}

function currentParameterRelationVersion() {
  return PARAMETER_RELATIONS.versions?.[state.versionId] || { programs: {} };
}

function parameterRelationsForNode(node) {
  if (!node || node.kind === "data") return { incompatible: [], requires: [], conditional: [], enables: [] };
  const programs = currentParameterRelationVersion().programs || {};
  const relations = programs[node.command] || programs[node.refId] || {};
  return {
    incompatible: relations.incompatible || [],
    requires: relations.requires || [],
    conditional: relations.conditional || [],
    enables: relations.enables || []
  };
}

function relationParameters(relation) {
  return (relation?.parameters || []).map(normalizeParameterName).filter(Boolean);
}

function namesIntersect(a, b) {
  return a.some((name) => b.has(name));
}

function activeInputParameterNames(node) {
  if (!node || node.kind === "data") return [];
  return Object.values(ensureInputValues(node))
    .filter((input) => {
      const connected = incomingEdgesForArtifact(node, input.artifact).length > 0;
      return connected || String(input.value || "").trim();
    })
    .map((input) => normalizeParameterName(input.flag || `--${input.artifact}`))
    .filter(Boolean);
}

function inputParameterNamesForArtifact(node, artifact) {
  const input = ensureInputValues(node)[artifact] || { artifact, flag: inferInputFlag(node, artifact) };
  const names = new Set([
    input.flag,
    `--${artifact}`,
    artifact
  ].map(normalizeParameterName).filter(Boolean));
  const artifactKey = normalizedDataKey(artifact);
  const artifactWords = significantArtifactWords(artifact);
  const exactParameters = [];
  const fuzzyParameters = [];
  (node.parameters || []).forEach((param) => {
    const candidates = parameterNameCandidates(param);
    const keys = (param.keys || []).map(normalizedDataKey);
    const candidateKeys = [...candidates.map(normalizedDataKey), ...keys];
    const exact = candidateKeys.some((name) => name === artifactKey);
    const wordMatch = artifactWords.length && candidateKeys.some((name) => artifactWords.every((word) => name.includes(word)));
    if (exact) exactParameters.push(param);
    else if (wordMatch) fuzzyParameters.push(param);
  });
  (exactParameters.length ? exactParameters : fuzzyParameters).forEach((param) => {
    parameterNameCandidates(param).forEach((name) => names.add(name));
  });
  return [...names];
}

function inputArtifactIsActive(node, artifact) {
  const input = ensureInputValues(node)[artifact];
  const connected = incomingEdgesForArtifact(node, artifact).length > 0;
  return connected || Boolean(String(input?.value || "").trim());
}

function inputArtifactCompatibilityState(node, artifact) {
  if (!node || node.kind === "data") return { status: "", reason: "" };
  const targetNames = inputParameterNamesForArtifact(node, artifact);
  const activeInputNames = new Set(activeInputParameterNames(node));
  const active = inputArtifactIsActive(node, artifact);
  if (active) {
    return {
      status: "active",
      reason: "This input is currently connected or has a value."
    };
  }
  const required = findRequiredRelation(node, targetNames, activeInputNames) || inputRequiredBySelectedRequiredSet(node, artifact, activeInputNames);
  if (required) {
    return {
      status: "required",
      reason: required.reason || relationNeedText(required.relation, targetNames, activeInputNames)
    };
  }
  const incompatible = findIncompatibleRelation(node, targetNames, activeInputNames);
  if (incompatible) {
    return {
      status: "incompatible",
      reason: relationLockText(incompatible, targetNames, activeInputNames)
    };
  }
  return { status: "", reason: "" };
}

function findRequiredRelation(node, targetNames, selectedNames) {
  const targetSet = new Set(targetNames);
  const selectedSet = selectedNames instanceof Set ? selectedNames : new Set(selectedNames);
  const targetParam = (node.parameters || []).find((param) => parameterNameCandidates(param).some((name) => targetSet.has(name)));
  const relations = [
    ...(parameterRelationsForNode(node).requires || []),
    ...(parameterRelationsForNode(node).conditional || [])
  ];
  return relations.find((relation) => {
    const params = relationParameters(relation);
    if (!params.some((name) => targetSet.has(name))) return false;
    if (!params.some((name) => selectedSet.has(name) && !targetSet.has(name))) return false;
    return Boolean(targetParam?.required) || params.length <= 2;
  }) || null;
}

function inputRequiredBySelectedRequiredSet(node, artifact, selectedNames) {
  const targetNames = new Set(inputParameterNamesForArtifact(node, artifact));
  const targetParam = (node.parameters || []).find((param) => parameterNameCandidates(param).some((name) => targetNames.has(name)));
  if (!targetParam?.required || !selectedNames.size) return null;
  const hasSelectedRequiredInput = (node.requires || []).some((otherArtifact) => {
    if (otherArtifact === artifact || !inputArtifactIsActive(node, otherArtifact)) return false;
    const otherNames = new Set(inputParameterNamesForArtifact(node, otherArtifact));
    const otherParam = (node.parameters || []).find((param) => parameterNameCandidates(param).some((name) => otherNames.has(name)));
    return Boolean(otherParam?.required);
  });
  return hasSelectedRequiredInput
    ? { reason: "Required with the selected input set." }
    : null;
}

function relationNeedText(relation, targetNames, selectedNames) {
  const others = relationOtherParameterNames(relation, targetNames, selectedNames);
  const otherText = others.length ? others.join(", ") : "the selected input";
  const evidence = relation?.evidence?.snippet ? ` ${compactText(relation.evidence.snippet, 150)}` : "";
  return `Needed with ${otherText}.${evidence}`;
}

function findIncompatibleRelation(node, targetNames, selectedNames) {
  const targetSet = new Set(targetNames);
  const selectedSet = selectedNames instanceof Set ? selectedNames : new Set(selectedNames);
  return (parameterRelationsForNode(node).incompatible || []).find((relation) => {
    const params = relationParameters(relation);
    if (!params.some((name) => targetSet.has(name))) return false;
    return params.some((name) => selectedSet.has(name) && !targetSet.has(name));
  }) || null;
}

function relationOtherParameterNames(relation, targetNames, selectedNames) {
  const targetSet = new Set(targetNames);
  const selectedSet = selectedNames instanceof Set ? selectedNames : new Set(selectedNames);
  return relationParameters(relation)
    .filter((name) => !targetSet.has(name))
    .filter((name) => selectedSet.size ? selectedSet.has(name) : true);
}

function relationLockText(relation, targetNames, selectedNames) {
  const others = relationOtherParameterNames(relation, targetNames, selectedNames);
  const otherText = others.length ? others.join(", ") : "another selected parameter";
  const evidence = relation?.evidence?.snippet ? ` ${compactText(relation.evidence.snippet, 150)}` : "";
  return `Disabled because it is incompatible with ${otherText}.${evidence}`;
}

function artifactPathPlaceholder(artifact, artifactType = "") {
  const key = normalizedDataKey(`${artifact} ${artifactType}`);
  if (key.includes("bam")) return "path/to/sample.bam";
  if (key.includes("fastq")) return "path/to/reads.fastq.gz";
  if (key.includes("fasta") || key.includes("genes-fasta") || key.includes("contigs")) return "path/to/contigs.fa";
  if (key.includes("json")) return "path/to/config.json";
  if (key.includes("txt") || key.includes("table")) return "path/to/table.txt";
  if (key.includes("db")) return `path/to/${normalizedDataKey(artifact) || "dataset"}.db`;
  if (key.includes("svg")) return "path/to/figure.svg";
  if (key.includes("vcf")) return "path/to/variants.vcf";
  return `path/to/${normalizedDataKey(artifact) || "dataset"}`;
}

function enforceParameterCompatibility(node, preferredParamName = "") {
  if (!node || node.kind === "data") return;
  const selected = new Set(activeInputParameterNames(node));
  const params = node.parameters || [];
  const preferred = normalizeParameterName(preferredParamName);
  const preferredIndex = preferred
    ? params.findIndex((param) => parameterNameCandidates(param).includes(preferred))
    : -1;
  const order = [
    ...(preferredIndex >= 0 ? [preferredIndex] : []),
    ...params.map((_, index) => index).filter((index) => index !== preferredIndex)
  ];

  order.forEach((index) => {
    const param = params[index];
    if (!param?.enabled) return;
    const targetNames = parameterNameCandidates(param);
    if (!targetNames.length) return;
    const conflict = findIncompatibleRelation(node, targetNames, selected);
    if (conflict && !param.required) {
      param.enabled = false;
      return;
    }
    targetNames.forEach((name) => selected.add(name));
  });
}

function parameterCompatibilityState(node) {
  const activeInputs = new Set(activeInputParameterNames(node));
  const selected = new Set(activeInputs);
  (node.parameters || []).forEach((param) => {
    if (param?.enabled) parameterNameCandidates(param).forEach((name) => selected.add(name));
  });

  const locks = new Map();
  (node.parameters || []).forEach((param, index) => {
    const targetNames = parameterNameCandidates(param);
    if (!targetNames.length) return;
    const inputConflict = findIncompatibleRelation(node, targetNames, activeInputs);
    const selectedConflict = findIncompatibleRelation(node, targetNames, selected);
    if (inputConflict && (param.required || !param.enabled)) {
      locks.set(index, {
        relation: inputConflict,
        reason: relationLockText(inputConflict, targetNames, activeInputs)
      });
      return;
    }
    if (!param.enabled && selectedConflict) {
      locks.set(index, {
        relation: selectedConflict,
        reason: relationLockText(selectedConflict, targetNames, selected)
      });
    }
  });
  return { activeInputs, selected, locks };
}

function ensureInputValues(node) {
  if (!node.inputValues) node.inputValues = {};
  (node.requires || []).forEach((artifact) => {
    if (!node.inputValues[artifact]) {
      node.inputValues[artifact] = {
        artifact,
        flag: inferInputFlag(node, artifact),
        value: ""
      };
    }
  });
  return node.inputValues;
}

function incomingEdgesForArtifact(node, artifact) {
  return activeEdges().filter((edge) => {
    if (edge.to !== node.id) return false;
    const from = nodeById(edge.from);
    return edge.label === artifact || (from?.kind === "data" && from.refId === artifact);
  });
}

function resolvedInputValue(node, input) {
  const explicitValue = String(input.value || "").trim();
  if (explicitValue) return explicitValue;
  const edge = incomingEdgesForArtifact(node, input.artifact)[0];
  if (!edge) return "";
  const from = nodeById(edge.from);
  if (from?.kind === "data") return from.value || from.refId || from.title || input.artifact;
  return edge.label || input.artifact;
}

function selectedCommandInputs(node) {
  return Object.values(ensureInputValues(node))
    .map((input) => ({
      ...input,
      resolvedValue: resolvedInputValue(node, input)
    }))
    .filter((input) => input.resolvedValue);
}

function inferInputFlag(node, artifact) {
  const artifactPattern = artifact.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
  for (const example of node.examples || []) {
    const match = String(example).match(new RegExp(`(^|\\s)(--[A-Za-z0-9-]+|-[A-Za-z])\\s+${artifactPattern}(\\s|$)`));
    if (match) return match[2];
  }
  const exact = (node.parameters || []).find((param) => param.name === `--${artifact}`);
  if (exact) return exact.name;
  const artifactWords = significantArtifactWords(artifact);
  const matchingParam = artifactWords.length
    ? (node.parameters || []).find((param) => {
      const candidates = [
        param.name,
        ...(param.aliases || []),
        ...(param.keys || [])
      ].map(normalizedDataKey);
      return candidates.some((name) => artifactWords.every((word) => name.includes(word)));
    })
    : null;
  if (matchingParam) return matchingParam.name;
  if (artifact.endsWith("-fasta") || artifact === "fasta") return "--fasta";
  if (artifact.endsWith("-db")) return `--${artifact}`;
  return `--${artifact}`;
}

function significantArtifactWords(artifact) {
  const generic = new Set(["db", "txt", "tsv", "csv", "json", "file", "data", "input", "output", "profile"]);
  return String(artifact || "")
    .toLowerCase()
    .split("-")
    .filter((word) => word.length > 2 && !generic.has(word));
}

function pushHistory() {
  state.history.past.push(JSON.stringify({
    tabs: state.tabs,
    activeTabId: state.activeTabId,
    settings: state.settings
  }));
  if (state.history.past.length > HISTORY_LIMIT) state.history.past.shift();
  state.history.future = [];
}

function restoreSnapshot(snapshot) {
  const parsed = JSON.parse(snapshot);
  state.tabs = parsed.tabs || [];
  state.activeTabId = parsed.activeTabId || state.tabs[0]?.id || "";
  state.settings = parsed.settings || state.settings;
  state.selectedNodeId = "";
  state.selectedEdgeId = "";
  state.focusedArtifact = null;
  state.connectStartId = "";
  applySettings();
  render();
}

function withHistory(mutator) {
  pushHistory();
  mutator();
  render();
}

function undo() {
  const snapshot = state.history.past.pop();
  if (!snapshot) return;
  state.history.future.push(JSON.stringify({
    tabs: state.tabs,
    activeTabId: state.activeTabId,
    settings: state.settings
  }));
  restoreSnapshot(snapshot);
}

function redo() {
  const snapshot = state.history.future.pop();
  if (!snapshot) return;
  state.history.past.push(JSON.stringify({
    tabs: state.tabs,
    activeTabId: state.activeTabId,
    settings: state.settings
  }));
  restoreSnapshot(snapshot);
}

function createTab(name = "Canvas", versionId = null) {
  return {
    id: uid("tab"),
    name,
    versionId: versionId || state.versionId,
    nodes: [],
    edges: []
  };
}

function init() {
  if (!DATABASE || !Array.isArray(DATABASE.versions)) {
    if (dom.sourceStatus) dom.sourceStatus.textContent = "Documentation database failed to load.";
    return;
  }
  state.db = applyLocalRuntimeDatabase(DATABASE, LOCAL_RUNTIME);

  const savedSettings = JSON.parse(localStorage.getItem("anvio-builder-settings") || "{}");
  state.settings = {
    ...state.settings,
    ...savedSettings
  };
  state.versionId = state.db.versions[0]?.id || "main";
  rebuildIntentGraphCache();
  state.tabs = [createTab("Canvas 1")];
  state.activeTabId = state.tabs[0].id;

  applySettings();
  bindEvents();
  render();
  requestAnimationFrame(centerCanvasViewport);
}

function applySettings() {
  document.documentElement.style.setProperty("--font-scale", state.settings.fontScale);
  document.body.classList.toggle("night", Boolean(state.settings.nightMode));
  applyPanelState();
  applyCanvasZoom();
  localStorage.setItem("anvio-builder-settings", JSON.stringify(state.settings));
}

function bindEvents() {
  document.addEventListener("click", handleActionClick);
  document.addEventListener("pointerdown", handleOutsideSelectionPointerDown);
  document.addEventListener("keydown", handleKeydown);
  document.querySelectorAll("[data-theme-switch]").forEach((switchInput) => {
    switchInput.addEventListener("change", () => {
      setNightMode(switchInput.checked);
    });
  });
  dom.searchInput.addEventListener("input", handleLibrarySearchInput);

  document.querySelectorAll(".segmented button").forEach((button) => {
    button.addEventListener("click", () => {
      document.querySelectorAll(".segmented button").forEach((other) => other.classList.remove("active"));
      button.classList.add("active");
      state.searchKind = button.dataset.kind;
      renderSearch();
    });
  });

  dom.resultList.addEventListener("click", (event) => {
    const card = event.target.closest("[data-spawn-kind]");
    if (!card) return;
    const item = findLibraryItem(card.dataset.spawnKind, card.dataset.spawnId);
    if (!item) return;
    state.focusedArtifact = null;
    const point = visibleSpawnPoint();
    withHistory(() => {
      spawnItem(item, point.x, point.y, true);
    });
  });

  dom.tabs.addEventListener("click", (event) => {
    if (event.target.closest("[data-tab-name-input]")) return;
    const tabButton = event.target.closest("[data-tab-id]");
    if (!tabButton) return;
    if (event.detail > 1) {
      event.preventDefault();
      beginTabNameEdit(tabButton.dataset.tabId);
      return;
    }
    if (state.activeTabId === tabButton.dataset.tabId) return;
    state.activeTabId = tabButton.dataset.tabId;
    state.selectedNodeId = "";
    state.selectedEdgeId = "";
    state.focusedArtifact = null;
    render();
  });
  dom.tabs.addEventListener("dblclick", handleTabDoubleClick);
  dom.tabs.addEventListener("keydown", handleTabEditKeydown);
  dom.tabs.addEventListener("focusout", handleTabEditFocusout);

  dom.nodeLayer.addEventListener("pointerdown", handleNodePointerDown);
  dom.nodeLayer.addEventListener("click", handleNodeLayerClick);
  dom.edgeLayer.addEventListener("pointerdown", handleEdgePointerDown);
  dom.canvasWorld.addEventListener("pointerdown", handleCanvasWorldPointerDown);
  dom.canvasFrame.addEventListener("pointerdown", handleCanvasPanPointerDown);
  dom.canvasFrame.addEventListener("wheel", handleCanvasWheel, { passive: false });
  dom.intentBuilderInput?.addEventListener("keydown", handleIntentBuilderKeydown);
  dom.intentBuilderInput?.addEventListener("input", handleIntentBuilderInput);

  dom.edgeLayer.addEventListener("pointermove", (event) => {
    const newHover = event.target.closest("[data-edge-id]")?.dataset.edgeId || null;

    if (newHover !== state.hoveredEdgeId) {
      state.hoveredEdgeId = newHover;
      state.mousePos = canvasPoint(event);
      renderEdges();
    } else if (state.hoveredEdgeId) {
      // Update mouse position while hovering, re-render controls
      state.mousePos = canvasPoint(event);
      renderEdges();
    }
  });

  dom.edgeLayer.addEventListener("pointerout", (event) => {
    if (event.relatedTarget && dom.edgeLayer.contains(event.relatedTarget)) return;
    if (state.hoveredEdgeId) {
      state.hoveredEdgeId = null;
      renderEdges();
    }
  });
  dom.edgeLayer.addEventListener("dblclick", (event) => {
    const edgeId = event.target.dataset.edgeId;
    if (!edgeId) return;
    event.preventDefault();
    withHistory(() => {
      const edge = activeEdges().find((candidate) => candidate.id === edgeId);
      if (edge) {
        edge.points = [];
        delete edge.layoutRoute;
      }
    });
  });

  dom.inspector.addEventListener("focusin", (event) => {
    if (!event.target.closest("[data-editable]") || textEditSnapshotOpen) return;
    pushHistory();
    textEditSnapshotOpen = true;
  });
  dom.inspector.addEventListener("focusout", () => {
    textEditSnapshotOpen = false;
  });
  dom.inspector.addEventListener("input", handleInspectorInput);
  dom.inspector.addEventListener("change", handleInspectorChange);

  dom.workspaceFileInput.addEventListener("change", handleWorkspaceFile);
  dom.workflowFileInput.addEventListener("change", handleWorkflowFile);
  dom.modal.addEventListener("click", handleModalBackdropClick);
  dom.modalContent.addEventListener("change", handleModalChange);
  dom.modalContent.addEventListener("input", handleModalInput);
}

function handleActionClick(event) {
  const button = event.target.closest("[data-action]");
  if (!button) return;
  const action = button.dataset.action;

  const actions = {
    "save-workspace": saveWorkspace,
    "load-workspace": () => dom.workspaceFileInput.click(),
    "import-workflow": openImportWorkflowDialog,
    "show-workflow-concepts": showWorkflowConcepts,
    "add-documented-pipeline": addDocumentedPipelineFromButton,
    "add-documented-pipeline-new-tab": addDocumentedPipelineToNewTabFromButton,
    "run-concept-search": runConceptSearchFromDialog,
    "select-workflow-concept": selectWorkflowConceptFromButton,
    undo,
    redo,
    "font-down": () => changeFont(-0.08),
    "font-up": () => changeFont(0.08),
    "zoom-out": () => zoomCanvas(-ZOOM_OUT_STEP),
    "zoom-in": () => zoomCanvas(ZOOM_IN_STEP),
    "center-canvas": centerCanvasViewport,
    "optimize-layout": optimizePlacement,
    "toggle-intent-builder": toggleIntentBuilder,
    "close-intent-builder": hideIntentBuilder,
    "run-intent-builder": runIntentBuilder,
    "toggle-library-panel": () => togglePanel("library"),
    "toggle-inspector-panel": () => togglePanel("inspector"),
    "toggle-night": toggleNightMode,
    help: showHelp,
    citation: showCitation,
    github: () => window.open(GITHUB_URL, "_blank", "noopener"),
    validate: showValidation,
    "new-tab": addTab,
    "remove-tab": () => removeTab(button.dataset.tabId),
    "export-open": openExportDialog,
    "export-format": exportSelectedFormat,
    "delete-node": deleteSelectedNode,
    "delete-edge": deleteSelectedEdge,
    "add-param": addCustomParam,
    "add-variable": addVariable,
    "reset-edge": resetSelectedEdge,
    "toggle-advanced": toggleAdvancedParams,
    "open-official-doc": openOfficialDocPopup,
    "confirm-import-workflow": confirmImportWorkflow
  };

  if (actions[action]) {
    event.preventDefault();
    actions[action](event);
  }
}

function handleKeydown(event) {
  const target = event.target;
  const isEditing = target && ["INPUT", "TEXTAREA", "SELECT"].includes(target.tagName);
  if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === "z") {
    event.preventDefault();
    if (event.shiftKey) redo();
    else undo();
  }
  if ((event.ctrlKey || event.metaKey) && event.key.toLowerCase() === "y") {
    event.preventDefault();
    redo();
  }
  if (!isEditing && (event.key === "Delete" || event.key === "Backspace")) {
    if (state.selectedNodeId) deleteSelectedNode();
    if (state.selectedEdgeId) deleteSelectedEdge();
  }
  if (event.key === "Escape") {
    if (dom.intentBuilder && !dom.intentBuilder.hidden) {
      hideIntentBuilder();
      if (isEditing) event.preventDefault();
    }
    state.mode = "select";
    state.connectStartId = "";
    renderToolbarState();
  }
}

function changeFont(delta) {
  state.settings.fontScale = Math.max(0.82, Math.min(1.28, Number((state.settings.fontScale + delta).toFixed(2))));
  applySettings();
}

function toggleNightMode() {
  setNightMode(!state.settings.nightMode);
}

function setNightMode(enabled) {
  state.settings.nightMode = Boolean(enabled);
  applySettings();
  renderToolbarState();
}

function revealInspectorPanel() {
  if (!state.settings.inspectorCollapsed) return;
  state.settings.inspectorCollapsed = false;
  applyPanelState();
  localStorage.setItem("anvio-builder-settings", JSON.stringify(state.settings));
  renderToolbarState();
}

function applyPanelState() {
  dom.workspace?.classList.toggle("library-collapsed", Boolean(state.settings.libraryCollapsed));
  dom.workspace?.classList.toggle("inspector-collapsed", Boolean(state.settings.inspectorCollapsed));
}

function togglePanel(which) {
  if (which === "library") {
    state.settings.libraryCollapsed = !state.settings.libraryCollapsed;
  }
  if (which === "inspector") {
    const hasSelection = Boolean(state.selectedNodeId || state.selectedEdgeId);
    state.settings.inspectorCollapsed = hasSelection ? !state.settings.inspectorCollapsed : true;
  }
  applySettings();
  renderToolbarState();
}

function toggleAdvancedParams(event) {
  const button = event.target.closest('[data-action="toggle-advanced"]');
  if (!button) return;
  const content = button.parentElement.querySelector('.advanced-content');
  if (!content) return;
  const isHidden = content.hasAttribute('hidden');
  if (isHidden) {
    content.removeAttribute('hidden');
    button.querySelector('.toggle-arrow').textContent = '\u25BC';
  } else {
    content.setAttribute('hidden', '');
    button.querySelector('.toggle-arrow').textContent = '\u25B6';
  }
}

function toggleConnectMode() {
  state.mode = state.mode === "connect" ? "select" : "connect";
  state.connectStartId = "";
  renderToolbarState();
}

function addTab() {
  withHistory(() => {
    const tab = createTab(`Canvas ${state.tabs.length + 1}`);
    state.tabs.push(tab);
    state.activeTabId = tab.id;
    state.selectedNodeId = "";
    state.selectedEdgeId = "";
  });
  requestAnimationFrame(centerCanvasViewport);
}

function removeTab(tabIdToRemove) {
  if (state.tabs.length < 2) {
    withHistory(() => {
      const tab = currentTab();
      if (!tab) return;
      tab.nodes = [];
      tab.edges = [];
      state.selectedNodeId = "";
      state.selectedEdgeId = "";
      state.focusedArtifact = null;
      state.connectStartId = "";
      state.mode = "select";
    });
    requestAnimationFrame(centerCanvasViewport);
    return;
  }
  withHistory(() => {
    let targetTabId = tabIdToRemove || state.activeTabId;
    const index = state.tabs.findIndex((tab) => tab.id === targetTabId);
    if (index === -1) return;
    state.tabs.splice(index, 1);
    if (state.activeTabId === targetTabId) {
      state.activeTabId = state.tabs[Math.max(0, index - 1)].id;
    }
    state.selectedNodeId = "";
    state.selectedEdgeId = "";
  });
}

function optimizePlacement() {
  if (!activeNodes().length) return;
  withHistory(() => {
    applyOptimizedPlacement(activeGraphCenter(viewportWorldCenter()));
  });
}

function toggleIntentBuilder() {
  if (!dom.intentBuilder) return;
  if (dom.intentBuilder.hidden) showIntentBuilder();
  else hideIntentBuilder();
}

function showIntentBuilder() {
  if (!dom.intentBuilder) return;
  dom.intentBuilder.hidden = false;
  updateIntentBuilderPreview("");
  requestAnimationFrame(() => dom.intentBuilderInput?.focus());
}

function hideIntentBuilder() {
  if (!dom.intentBuilder) return;
  dom.intentBuilder.hidden = true;
}

function handleIntentBuilderKeydown(event) {
  if (event.key === "Enter") {
    event.preventDefault();
    runIntentBuilder();
  }
  if (event.key === "Escape") {
    event.preventDefault();
    hideIntentBuilder();
  }
}

function handleIntentBuilderInput() {
  clearTimeout(intentPreviewTimer);
  if (dom.intentBuilderPreview) {
    dom.intentBuilderPreview.textContent = "";
  }
}

function updateIntentBuilderPreview(html) {
  if (!dom.intentBuilderPreview) return;
  dom.intentBuilderPreview.innerHTML = html;
}

function runIntentBuilder() {
  const query = String(dom.intentBuilderInput?.value || "").trim();
  if (!query) {
    updateIntentBuilderPreview("Describe a goal first, for example: <strong>profile metagenomic reads against a contigs database</strong>.");
    return;
  }

  const graph = cachedIntentGraph();
  const context = intentQueryContext(query, graph);
  const plan = planGraphIntentWorkflowFromKeywords(context, graph);
  if (plan.error) {
    updateIntentBuilderPreview(plan.error);
    return;
  }

  let result = null;
  withHistory(() => {
    result = buildGraphIntentWorkflow(plan, context, graph);
  });

  const inputText = result.createdInputs.length
    ? `Added input data: ${result.createdInputs.map(escapeHtml).join(", ")}.`
    : "No new input data block was added.";
  const programText = result.createdPrograms.length
    ? ` Added program blocks: ${result.createdPrograms.map(escapeHtml).join(", ")}.`
    : "";
  const reusedText = result.reusedInputs.length
    ? ` Reused existing input blocks: ${result.reusedInputs.map(escapeHtml).join(", ")}.`
    : "";
  const ambiguityText = result.ambiguousInputs.length
    ? ` Ambiguous inputs left for you to choose: ${result.ambiguousInputs.map(escapeHtml).join(", ")}.`
    : "";
  updateIntentBuilderPreview(`<strong>${escapeHtml(result.title)}</strong> built from ${escapeHtml(graph.source || "the cached graph")}.${programText} ${inputText}${reusedText}${ambiguityText}`);
}

function previewIntentGraphPlan(query) {
  if (!query) {
    updateIntentBuilderPreview("");
    return;
  }
  const goal = inferIntentGoal(query);
  if (!goal) {
    updateIntentBuilderPreview("No internal anvi'o graph endpoint matched yet.");
    return;
  }
  const context = intentQueryContext(query);
  const inputNames = context.mentionedArtifacts.slice(0, 4).map((artifact) => artifact.id);
  const pathText = [
    ...inputNames.map((name) => `<span class="chip">${escapeHtml(name)}</span>`),
    `<span class="chip relation">${escapeHtml(intentGoalTitle(goal))}</span>`
  ].join(" ");
  updateIntentBuilderPreview(`
    <strong>${escapeHtml(intentGoalTitle(goal))}</strong> is the current endpoint.
    ${pathText ? `<span class="intent-path">${pathText}</span>` : ""}
    <span class="intent-hint">Press Build to compute the optimal path.</span>
  `);
}

function inferIntentGoal(query) {
  const graph = cachedIntentGraph();
  const context = intentQueryContext(query, graph);

  const programCandidates = graph.programSearch
    .map((doc) => {
      const program = doc.program;
      let score = context.queryMagnitude ? simpleVectorCosine(context.queryVector, doc.vector, context.queryMagnitude, doc.magnitude) * 500 : 0;
      score += scoreField(program.command, context.normalized, context.words, {
        exact: 1200,
        prefix: 720,
        phrase: 420,
        wordExact: 90,
        wordPartial: 48
      });
      score += scoreField(program.description, context.normalized, context.words, {
        exact: 120,
        prefix: 70,
        phrase: 45,
        wordExact: 16,
        wordPartial: 8
      });
      context.mentionedArtifacts.forEach((artifact) => {
        if ((program.requires || []).some((required) => artifactCompatible(required, artifact.id))) score += 180;
        if ((program.provides || []).some((provided) => artifactCompatible(provided, artifact.id))) score += 140;
      });
      return { type: "program", program, score };
    })
    .filter((candidate) => candidate.score > 0)
    .sort((a, b) => b.score - a.score || (a.program.command || "").localeCompare(b.program.command || ""));

  const artifactCandidates = graph.artifactSearch
    .map((doc) => ({
      type: "artifact",
      artifact: doc.artifact,
      score: scoreIntentArtifact(doc, context)
    }))
    .filter((candidate) => candidate.score > 0 && (graph.producersByArtifact.get(candidate.artifact.id)?.length || context.mentionedIds.has(candidate.artifact.id)))
    .sort((a, b) => b.score - a.score || a.artifact.id.localeCompare(b.artifact.id));

  const bestProgram = programCandidates[0] || null;
  const bestArtifact = artifactCandidates[0] || null;
  if (bestArtifact && (!bestProgram || bestArtifact.score > bestProgram.score * 1.18)) return bestArtifact;
  return bestProgram;
}

function simpleVectorCosine(queryVector, docVector, queryMagnitude, docMagnitude) {
  if (!queryMagnitude || !docMagnitude) return 0;
  let dot = 0;
  queryVector.forEach((value, token) => {
    dot += value * (docVector.get(token) || 0);
  });
  return dot / (queryMagnitude * docMagnitude);
}

function scoreIntentArtifact(doc, context) {
  const artifact = doc.artifact;
  let score = context.queryMagnitude ? simpleVectorCosine(context.queryVector, doc.vector, context.queryMagnitude, doc.magnitude) * 420 : 0;
  score += scoreField(artifact.id, context.normalized, context.words, {
    exact: 1600,
    prefix: 900,
    phrase: 520,
    wordExact: 90,
    wordPartial: 44
  });
  score += scoreField(artifact.title, context.normalized, context.words, {
    exact: 1300,
    prefix: 720,
    phrase: 420,
    wordExact: 70,
    wordPartial: 34
  });
  if (/\bmake|create|generate|build|produce|obtenir|generer|creer|faire\b/i.test(context.query)) score += 80;
  return score;
}

function intentQueryContext(query, graph = cachedIntentGraph()) {
  const normalized = normalizedIntentText(query);
  const segments = intentQuerySegments(query);
  const queryVector = textVector(query);
  const words = tokenizeConceptText(query);
  const sourceMatch = segmentMatchContext(segments.source || query);
  const goalMatch = segmentMatchContext(segments.goal || query);
  const sourceCandidates = sourceArtifactCandidates(segments.source || query, graph);
  const explicitInputArtifacts = sourceCandidates.map((candidate) => candidate.artifact);
  const mentioned = mergeArtifactMentions(mentionedArtifacts(query, graph, normalized), explicitInputArtifacts);
  return {
    query,
    normalized,
    segments,
    words,
    uniqueWords: [...new Set(words)],
    queryVector,
    queryMagnitude: vectorMagnitude(queryVector),
    sourceMatch,
    goalMatch,
    mentionedArtifacts: mentioned,
    mentionedIds: new Set(mentioned.map((artifact) => artifact.id)),
    sourceArtifactScores: new Map(sourceCandidates.map((candidate) => [candidate.artifact.id, candidate.score])),
    explicitInputArtifacts,
    explicitInputIds: new Set(explicitInputArtifacts.map((artifact) => artifact.id))
  };
}

function mentionedArtifacts(query, graph = cachedIntentGraph(), normalized = normalizedIntentText(query)) {
  return graph.artifactSearch
    .filter((doc) => doc.mentionKeys.some((candidate) => normalized.includes(candidate)))
    .map((doc) => doc.artifact)
    .sort((a, b) => a.id.localeCompare(b.id));
}

function intentQuerySegments(query) {
  const normalized = normalizedIntentText(query);
  const goalMatch = normalized.match(/\b(i want to|i need to|want to|need to|run|perform|create|generate|build|produce|obtenir|generer|creer|faire|lancer)\b/);
  if (!goalMatch) return { source: normalized, goal: normalized };
  return {
    source: normalized.slice(0, goalMatch.index).trim(),
    goal: normalized.slice(goalMatch.index).trim()
  };
}

function sourceArtifactCandidates(sourceText, graph = cachedIntentGraph()) {
  const context = segmentMatchContext(sourceText);
  if (!context.uniqueWords.length) return [];
  const scored = graph.artifactSearch
    .map((doc) => ({
      artifact: doc.artifact,
      score: scoreArtifactAgainstSegment(doc, context)
    }))
    .filter((candidate) => candidate.score > 0)
    .sort((a, b) => b.score - a.score || a.artifact.id.localeCompare(b.artifact.id));
  if (!scored.length) return [];
  const best = scored[0].score;
  return scored
    .filter((candidate) => candidate.score >= Math.max(90, best * 0.82))
    .slice(0, 6);
}

function segmentMatchContext(text) {
  const normalized = normalizedIntentText(text);
  const vector = textVector(text);
  return {
    query: text,
    normalized,
    words: tokenizeConceptText(text),
    uniqueWords: [...new Set(tokenizeConceptText(text))],
    queryVector: vector,
    queryMagnitude: vectorMagnitude(vector)
  };
}

function scoreArtifactAgainstSegment(doc, context) {
  const artifact = doc.artifact;
  const keywordWeights = doc.keywordWeights || new Map();
  let score = 0;
  let matched = 0;
  context.uniqueWords.forEach((word) => {
    const weight = keywordWeights.get(word) || 0;
    if (weight) {
      score += weight * 20;
      matched += 1;
    } else if (doc.searchText?.includes(word)) {
      score += 10;
    }
  });
  score += simpleVectorCosine(context.queryVector, doc.vector, context.queryMagnitude, doc.magnitude) * 260;
  score += scoreField(artifactKeywordText(artifact), context.normalized, context.uniqueWords, {
    exact: 1200,
    prefix: 680,
    phrase: 420,
    wordExact: 110,
    wordPartial: 48
  });
  if (!matched) score *= 0.35;
  else score *= 0.8 + matched / Math.max(1, context.uniqueWords.length);
  return score;
}

function mergeArtifactMentions(...groups) {
  const byId = new Map();
  groups.flat().forEach((artifact) => {
    if (!artifact?.id) return;
    if (!byId.has(artifact.id)) byId.set(artifact.id, artifact);
  });
  return [...byId.values()].sort((a, b) => a.id.localeCompare(b.id));
}

function normalizedIntentText(value) {
  return String(value || "")
    .toLowerCase()
    .normalize("NFKD")
    .replace(/[\u0300-\u036f]/g, "")
    .replace(/[^a-z0-9]+/g, " ")
    .trim();
}

function buildGraphIntentWorkflow(plan, queryOrContext, graph = cachedIntentGraph()) {
  const tab = currentTab();
  const center = visibleSpawnPoint();
  const context = typeof queryOrContext === "string" ? intentQueryContext(queryOrContext, graph) : queryOrContext;
  const nodesByProgram = new Map();
  const existingProgramIds = new Set(activeNodes().filter((node) => node.kind === "program").map(programNodeKey));
  const createdPrograms = [];
  const createdInputs = [];
  const reusedInputs = [];
  const ambiguousInputs = [];

  plan.programs.forEach((plannedProgram, index) => {
    const key = programItemKey(plannedProgram);
    const existing = activeNodes().find((node) => node.kind === "program" && programNodeKey(node) === key);
    const node = existing || createNodeFromItem(plannedProgram, center.x + index * 320, center.y);
    if (!existing) {
      tab.nodes.push(node);
      createdPrograms.push(node.command || node.title);
    }
    nodesByProgram.set(key, node);
  });

  plan.programs.forEach((plannedProgram) => {
    const consumer = nodesByProgram.get(programItemKey(plannedProgram));
    if (!consumer) return;
    graphInputArtifacts(plannedProgram, context, graph).forEach((artifactId) => {
      const plannedProducer = findPlannedProducerNode(artifactId, consumer.id, nodesByProgram);
      if (plannedProducer) {
        ensureEdge(plannedProducer.id, consumer.id, artifactId);
        return;
      }

      const existingProducer = existingProducerForArtifact(artifactId, consumer.id);
      if (existingProducer) {
        ensureEdge(existingProducer.id, consumer.id, artifactId);
        reusedInputs.push(artifactId);
        return;
      }

      const dataItem = findDataBlock(artifactId);
      const dataNode = reusableDataNode(dataItem, consumer.x - 330, consumer.y);
      if (dataNode && !dataNodeHasProducer(dataNode)) {
        ensureEdge(dataNode.id, consumer.id, artifactId);
        reusedInputs.push(artifactId);
        return;
      }

      if (!plan.initialArtifacts.has(artifactId)) {
        ambiguousInputs.push(artifactId);
        return;
      }

      const y = laneY(consumer.y, createdInputs.length, Math.max(1, plan.initialArtifacts.size));
      const newDataNode = createNodeFromItem(dataItem, consumer.x - 330, y);
      newDataNode.value = artifactId;
      tab.nodes.push(newDataNode);
      ensureEdge(newDataNode.id, consumer.id, artifactId);
      createdInputs.push(artifactId);
    });
  });

  plan.programs.forEach((plannedProgram) => {
    connectExistingFunctionBlocks(nodesByProgram.get(programItemKey(plannedProgram)));
  });

  const node = plan.goal?.type === "program"
    ? nodesByProgram.get(programItemKey(plan.goal.program)) || activeNodes().find((candidate) => candidate.kind === "program" && candidate.command === plan.goal.program.command)
    : nodesByProgram.get(programItemKey(plan.programs[plan.programs.length - 1]));
  state.selectedNodeId = node?.id || "";
  state.selectedEdgeId = "";
  state.focusedArtifact = null;
  revealInspectorPanel();
  applyOptimizedPlacement(activeGraphCenter(center));
  return {
    node,
    title: plan.title || intentGoalTitle(plan.goal),
    createdPrograms: createdPrograms.filter((name) => !existingProgramIds.has(name)),
    createdInputs: [...new Set(createdInputs)],
    reusedInputs: [...new Set(reusedInputs)],
    ambiguousInputs: [...new Set(ambiguousInputs)]
  };
}

function planGraphIntentWorkflow(goal, queryOrContext) {
  const graph = cachedIntentGraph();
  const context = typeof queryOrContext === "string" ? intentQueryContext(queryOrContext, graph) : queryOrContext;
  const programs = new Map();
  const initialArtifacts = new Set();
  const visiting = new Set();
  const plannedOrder = [];
  const mentionedEntryArtifacts = context.mentionedIds;

  function addProgram(program, depth = 0) {
    const key = programItemKey(program);
    if (!key || visiting.has(key)) return;
    if (programs.has(key)) return;
    visiting.add(key);

    graphInputArtifacts(program, context, graph).forEach((artifactId) => {
      if (artifactSatisfiedByExistingOrPlanned(artifactId, programs)) return;
      if (hasCompatibleArtifactInSet(artifactId, mentionedEntryArtifacts)) {
        initialArtifacts.add(artifactId);
        return;
      }
      const producer = depth < 4 ? bestProducerProgramForArtifact(artifactId, context, programs, graph) : null;
      if (producer) addProgram(producer, depth + 1);
      if (!artifactSatisfiedByExistingOrPlanned(artifactId, programs) && !context.explicitInputIds?.size) initialArtifacts.add(artifactId);
    });

    visiting.delete(key);
    programs.set(key, program);
    plannedOrder.push(program);
  }

  if (goal.type === "program") {
    addProgram(goal.program, 0);
  } else {
    const producer = bestProducerProgramForArtifact(goal.artifact.id, context, programs, graph);
    if (producer) addProgram(producer, 0);
    else initialArtifacts.add(goal.artifact.id);
  }
  return { programs: plannedOrder, initialArtifacts };
}

function artifactSatisfiedByExistingOrPlanned(artifactId, plannedPrograms) {
  if (existingProducerForArtifact(artifactId)) return true;
  if (activeNodes().some((node) => node.kind === "data" && artifactCompatible(artifactId, node.refId))) return true;
  return [...plannedPrograms.values()].some((program) => programProvidesArtifact(program, artifactId));
}

function findPlannedProducerNode(artifactId, consumerId, nodesByProgram) {
  return [...nodesByProgram.values()]
    .filter((node) => node.id !== consumerId && node.kind === "program" && nodeProvidesArtifact(node, artifactId))
    .sort(compareLayoutNodes)[0] || null;
}

function bestProducerProgramForArtifact(artifactId, queryOrContext, plannedPrograms = new Map(), graph = cachedIntentGraph()) {
  const context = typeof queryOrContext === "string" ? intentQueryContext(queryOrContext, graph) : queryOrContext;
  const candidates = (graph.producersByArtifact.get(artifactId) || [])
    .filter((program) => !plannedPrograms.has(programItemKey(program)))
    .map((program) => ({
      program,
      score: scoreIntentProducer(program, artifactId, context, plannedPrograms, graph)
    }))
    .filter((candidate) => candidate.score > 0)
    .sort((a, b) => b.score - a.score || graphInputArtifacts(a.program, context, graph).length - graphInputArtifacts(b.program, context, graph).length);

  if (!candidates.length) return null;
  if (candidates.length === 1) return candidates[0].program;
  const [best, second] = candidates;
  return best.score >= second.score + 12 ? best.program : null;
}

function scoreIntentProducer(program, artifactId, queryOrContext, plannedPrograms = new Map(), graph = cachedIntentGraph()) {
  const context = typeof queryOrContext === "string" ? intentQueryContext(queryOrContext, graph) : queryOrContext;
  const exactOutput = (program.provides || []).some((provided) => artifactCompatible(artifactId, provided));
  const doc = graph.programSearchByKey.get(programItemKey(program));
  const similarity = doc ? simpleVectorCosine(context.queryVector, doc.vector, context.queryMagnitude, doc.magnitude) : 0;
  const estimatedCost = estimateProgramGraphCost(program, context, plannedPrograms, graph);
  let score = exactOutput ? 140 : 90;
  score += similarity * 80;
  score += scoreField(program.command, context.normalized, context.words, {
    exact: 90,
    prefix: 55,
    phrase: 34,
    wordExact: 8,
    wordPartial: 4
  });
  score -= estimatedCost * 12;
  if (/^anvi-gen-/.test(program.command || "") && artifactId.endsWith("-db")) score += 8;
  if ((program.command || "").includes("merge") && !/\bmerge|merged|multi|multiple\b/i.test(context.query)) score -= 16;
  if ((program.command || "").includes("profile") && /\bread|reads|bam|mapping|metagenom/i.test(context.query)) score += 18;
  return score;
}

function graphInputArtifacts(program, queryOrContext = "", graph = cachedIntentGraph()) {
  const context = typeof queryOrContext === "string" ? intentQueryContext(queryOrContext, graph) : queryOrContext;
  const artifacts = [...new Set(program?.requires || [])].filter(Boolean);
  if (!artifacts.length) return [];
  const cached = graph.essentialInputsByProgram.get(programItemKey(program));
  const essential = cached ? [...cached] : baseGraphInputArtifacts(program);
  const mentioned = context.mentionedIds;
  artifacts.forEach((artifact) => {
    if (hasCompatibleArtifactInSet(artifact, mentioned) && !essential.includes(artifact)) essential.push(artifact);
  });
  return essential.sort((a, b) => inputArtifactPriority(a, context) - inputArtifactPriority(b, context) || a.localeCompare(b));
}

function inputArtifactPriority(artifactId, context) {
  const sourceScore = compatibleSourceArtifactScore(artifactId, context);
  if (sourceScore > 0) return 0;
  if (hasCompatibleArtifactInSet(artifactId, context.mentionedIds || new Set())) return 1;
  return 3;
}

function compatibleSourceArtifactScore(artifactId, context) {
  let best = 0;
  (context.sourceArtifactScores || new Map()).forEach((score, candidate) => {
    if (artifactCompatible(artifactId, candidate)) best = Math.max(best, score);
  });
  return best;
}

function baseGraphInputArtifacts(program) {
  const artifacts = [...new Set(program?.requires || [])].filter(Boolean);
  if (!artifacts.length) return [];
  const requiredParams = (program.parameters || []).filter((parameter) => parameter.required);
  if (!requiredParams.length) return artifacts;
  return artifacts.filter((artifact) => {
    const artifactKey = normalizedDataKey(artifact);
    return requiredParams.some((parameter) => parameterMatchesArtifact(parameter, artifactKey));
  });
}

function parameterMatchesArtifact(parameter, artifactKey) {
  const names = [
    parameter.name,
    ...(parameter.aliases || []),
    ...(parameter.keys || [])
  ].map(normalizedDataKey).filter(Boolean);
  return names.some((name) => name === artifactKey || name.includes(artifactKey) || artifactKey.includes(name));
}

function programItemKey(program) {
  return program?.command || program?.id || "";
}

function programNodeKey(node) {
  return node?.command || node?.refId || node?.title || "";
}

function programProvidesArtifact(program, artifactId) {
  return (program?.provides || []).some((provided) => artifactCompatible(artifactId, provided));
}

function planGraphIntentWorkflowFromKeywords(context, graph = cachedIntentGraph()) {
  const matches = graphKeywordMatches(context, graph);
  if (!matches.length) {
    return { error: intentQueryGuidance(context, "low") };
  }

  const explored = exploreIntentGraphNeighborhood(matches, graph, 3);
  const pathPlan = bestIntentPathPlan(matches, explored, context, graph);
  if (pathPlan) return pathPlan;

  const goal = inferIntentGoalFromMatches(matches, context, graph) || inferIntentGoal(context.query);
  if (goal) {
    const fallbackPlan = planGraphIntentWorkflow(goal, context);
    if (fallbackPlan.programs.length) {
      return {
        ...fallbackPlan,
        goal,
        title: intentGoalTitle(goal)
      };
    }
  }

  return { error: intentQueryGuidance(context, context.uniqueWords.length > 10 || matches.length > 35 ? "high" : "low") };
}

function graphKeywordMatches(context, graph) {
  const candidates = [
    ...graph.programSearch.map((doc) => ({ type: "program", nodeId: doc.nodeId, doc, item: doc.program })),
    ...graph.artifactSearch.map((doc) => ({ type: "artifact", nodeId: doc.nodeId, doc, item: doc.artifact }))
  ];

  return candidates
    .map((candidate) => {
      const sourceScore = graphKeywordRoleScore(candidate, context.sourceMatch || context, context, "source");
      const goalScore = graphKeywordRoleScore(candidate, context.goalMatch || context, context, "goal");
      return {
        ...candidate,
        sourceScore,
        goalScore,
        score: Math.max(sourceScore, goalScore)
      };
    })
    .filter((candidate) => candidate.score > 0)
    .sort((a, b) => b.score - a.score || intentMatchTitle(a).localeCompare(intentMatchTitle(b)))
    .slice(0, 48);
}

function graphKeywordScore(candidate, context) {
  return graphKeywordRoleScore(
    candidate,
    candidate.type === "program" ? (context.goalMatch || context) : (context.sourceMatch || context),
    context,
    candidate.type === "program" ? "goal" : "source"
  );
}

function graphKeywordRoleScore(candidate, matchContext, context, role) {
  const doc = candidate.doc;
  const item = candidate.item;
  const keywordSet = doc.keywordSet || new Set(doc.keywords || []);
  const keywordWeights = doc.keywordWeights || new Map();
  let score = 0;
  let matchedWords = 0;
  matchContext.uniqueWords.forEach((word) => {
    const weight = keywordWeights.get(word) || 0;
    if (weight) {
      score += 18 * weight;
      matchedWords += 1;
    } else if (keywordSet.has(word)) {
      score += 45;
      matchedWords += 1;
    } else if (doc.searchText?.includes(word)) {
      score += 12;
    }
  });
  if (matchContext.normalized && doc.searchText?.includes(matchContext.normalized)) score += 220;
  score += simpleVectorCosine(matchContext.queryVector, doc.vector, matchContext.queryMagnitude, doc.magnitude) * 260;
  if (matchContext.uniqueWords.length) {
    score *= 0.7 + (matchedWords / matchContext.uniqueWords.length) * 0.45;
  }

  if (candidate.type === "program") {
    score += scoreField(item.command, matchContext.normalized, matchContext.uniqueWords, {
      exact: 1500,
      prefix: 760,
      phrase: 420,
      wordExact: 120,
      wordPartial: 55
    });
  } else {
    const artifactText = artifactKeywordText(item).toLowerCase();
    score += scoreField(item.id, matchContext.normalized, matchContext.uniqueWords, {
      exact: 1500,
      prefix: 760,
      phrase: 420,
      wordExact: 120,
      wordPartial: 55
    });
    score += scoreField(artifactText, matchContext.normalized, matchContext.uniqueWords, {
      exact: 1000,
      prefix: 580,
      phrase: 340,
      wordExact: 90,
      wordPartial: 38
    });
    if (role === "source") score += (context.sourceArtifactScores?.get(item.id) || 0) * 1.4;
    if (role === "goal" && context.explicitInputIds?.has(item.id)) score *= 0.2;
  }
  return score;
}

function intentMatchTitle(match) {
  if (match.type === "program") return match.item.command || match.item.id || "";
  return match.item.id || match.item.title || "";
}

function exploreIntentGraphNeighborhood(matches, graph, maxDepth = 3) {
  const nodes = new Set();
  const edgeKeys = new Set();
  const queue = matches.map((match) => ({ nodeId: match.nodeId, depth: 0 }));
  matches.forEach((match) => nodes.add(match.nodeId));

  while (queue.length) {
    const current = queue.shift();
    if (current.depth >= maxDepth) continue;
    const edges = [
      ...(graph.outgoingByNode.get(current.nodeId) || []),
      ...(graph.incomingByNode.get(current.nodeId) || [])
    ];
    edges.forEach((edge) => {
      const next = edge.from === current.nodeId ? edge.to : edge.from;
      edgeKeys.add(`${edge.from}>${edge.to}:${edge.artifact || ""}`);
      if (nodes.has(next)) return;
      nodes.add(next);
      queue.push({ nodeId: next, depth: current.depth + 1 });
    });
  }
  return { nodes, edgeKeys };
}

function bestIntentPathPlan(matches, explored, context, graph) {
  const topMatches = matches.slice(0, 18);
  const sourceMatches = topMatches
    .filter((match) => match.sourceScore > 0)
    .sort((a, b) => b.sourceScore - a.sourceScore)
    .slice(0, 8);
  const targetMatches = topMatches
    .filter((match) => match.goalScore > 0 && !context.explicitInputIds?.has(match.item?.id))
    .sort((a, b) => b.goalScore - a.goalScore)
    .slice(0, 10);
  let best = null;
  sourceMatches.forEach((source) => {
    targetMatches.forEach((target) => {
      if (source.nodeId === target.nodeId) return;
      const path = directedShortestIntentPath(source.nodeId, target.nodeId, graph, explored.nodes, 8);
      if (!path || path.length < 2) return;
      const programs = programsFromGraphPath(path, graph);
      if (!programs.length) return;
      const coveredScore = topMatches
        .filter((match) => path.includes(match.nodeId))
        .filter((match) => match.nodeId !== source.nodeId && match.nodeId !== target.nodeId)
        .reduce((sum, match) => sum + Math.min(match.score, 260) * 0.12, 0);
      const score = coveredScore + source.sourceScore + target.goalScore * 1.35 - path.length * 42 - programs.length * 44;
      if (!best || score > best.score) {
        best = { path, programs, score, source, target };
      }
    });
  });
  if (!best) return null;
  return planFromIntentPath(best, context, graph);
}

function directedShortestIntentPath(sourceId, targetId, graph, allowedNodes, maxHops = 8) {
  const queue = [{ nodeId: sourceId, path: [sourceId] }];
  const seen = new Set([sourceId]);
  while (queue.length) {
    const current = queue.shift();
    if (current.path.length - 1 >= maxHops) continue;
    const edges = graph.outgoingByNode.get(current.nodeId) || [];
    for (const edge of edges) {
      if (!allowedNodes.has(edge.to)) continue;
      if (seen.has(edge.to)) continue;
      const path = [...current.path, edge.to];
      if (edge.to === targetId) return path;
      seen.add(edge.to);
      queue.push({ nodeId: edge.to, path });
    }
  }
  return null;
}

function programsFromGraphPath(path, graph) {
  const programs = [];
  const seen = new Set();
  path.forEach((nodeId) => {
    if (!nodeId.startsWith("program:")) return;
    const program = graph.nodes.get(nodeId)?.item;
    const key = programItemKey(program);
    if (!program || seen.has(key)) return;
    seen.add(key);
    programs.push(program);
  });
  return programs;
}

function planFromIntentPath(pathMatch, context, graph) {
  const programs = pathMatch.programs;
  const initialArtifacts = new Set();
  const plannedPrograms = new Map(programs.map((program) => [programItemKey(program), program]));
  const pathArtifacts = new Set(pathMatch.path.filter((nodeId) => nodeId.startsWith("artifact:")).map((nodeId) => nodeId.slice("artifact:".length)));

  programs.forEach((program) => {
    graphInputArtifacts(program, context, graph).forEach((artifactId) => {
      if (artifactSatisfiedByExistingOrPlanned(artifactId, plannedPrograms)) return;
      const explicitlyMentioned = hasCompatibleArtifactInSet(artifactId, context.mentionedIds);
      const sourceScore = compatibleSourceArtifactScore(artifactId, context);
      const pathOnlyInput = hasCompatibleArtifactInSet(artifactId, pathArtifacts);
      if (sourceScore > 0 || explicitlyMentioned || (pathOnlyInput && !context.explicitInputIds?.size)) {
        initialArtifacts.add(artifactId);
        return;
      }
      if (context.explicitInputIds?.size) {
        return;
      }
      if (!(graph.producersByArtifact.get(artifactId) || []).length) initialArtifacts.add(artifactId);
    });
  });

  return {
    programs,
    initialArtifacts,
    goal: graphGoalFromNodeId(pathMatch.target.nodeId, graph),
    title: intentPathTitle(pathMatch, graph),
    path: pathMatch.path
  };
}

function graphGoalFromNodeId(nodeId, graph) {
  const node = graph.nodes.get(nodeId);
  if (!node) return null;
  if (node.kind === "program") return { type: "program", program: node.item };
  return { type: "artifact", artifact: node.item };
}

function inferIntentGoalFromMatches(matches, context, graph) {
  const best = matches
    .map((match) => {
      let score = match.score;
      if (match.type === "artifact" && (graph.producersByArtifact.get(match.item.id) || []).length) score += 90;
      if (match.type === "program") score += 60;
      if (/\bmake|create|generate|build|produce|obtenir|generer|creer|faire\b/i.test(context.query) && match.type === "artifact") score += 80;
      return { match, score };
    })
    .sort((a, b) => b.score - a.score)[0]?.match;
  return best ? graphGoalFromNodeId(best.nodeId, graph) : null;
}

function intentPathTitle(pathMatch, graph) {
  const source = graph.nodes.get(pathMatch.source.nodeId)?.item;
  const target = graph.nodes.get(pathMatch.target.nodeId)?.item;
  const sourceTitle = source?.command || source?.id || source?.title || "graph match";
  const targetTitle = target?.command || target?.id || target?.title || "target";
  return `${sourceTitle} to ${targetTitle}`;
}

function intentQueryGuidance(context, density) {
  if (density === "high") {
    return "I found too many possible graph anchors but no coherent path. Try a shorter request with the starting data and final artifact or command.";
  }
  return "I could not find enough graph anchors to build a workflow. Please be more specific, for example by naming the input data and the expected final artifact.";
}

function cachedIntentGraph() {
  const version = currentVersion();
  const cacheKey = `${version.id}:${version.programs?.length || 0}:${version.dataBlocks?.length || 0}`;
  if (intentGraphCache.has(cacheKey)) return intentGraphCache.get(cacheKey);
  const graph = buildIntentGraph(version);
  intentGraphCache.clear();
  intentGraphCache.set(cacheKey, graph);
  return graph;
}

function rebuildIntentGraphCache() {
  const version = currentVersion();
  if (!version) return null;
  const cacheKey = `${version.id}:${version.programs?.length || 0}:${version.dataBlocks?.length || 0}`;
  const graph = buildIntentGraph(version);
  intentGraphCache.clear();
  intentGraphCache.set(cacheKey, graph);
  return graph;
}

function buildIntentGraph(version) {
  if (version.programNetwork?.nodes?.length && version.programNetwork?.links?.length) {
    return buildIntentGraphFromAnvioNetwork(version, version.programNetwork);
  }
  return buildIntentGraphFromRuntimeMetadata(version);
}

function buildIntentGraphFromAnvioNetwork(version, network) {
  const artifacts = new Map((version.dataBlocks || []).map((artifact) => [artifact.id, artifact]));
  const programs = new Map(
    (version.programs || [])
      .filter(isInternalProgramItem)
      .map((program) => [programItemKey(program), program])
  );
  programs.forEach((program) => {
    [...(program.requires || []), ...(program.provides || [])].forEach((artifactId) => {
      if (artifactId && !artifacts.has(artifactId)) artifacts.set(artifactId, findDataBlock(artifactId));
    });
  });
  (network.nodes || []).forEach((node) => {
    if (node.type === "PROGRAM") return;
    const id = String(node.id || node.name || "").toLowerCase();
    if (id && !artifacts.has(id)) {
      artifacts.set(id, {
        id,
        kind: "data",
        title: node.name || id,
        artifactType: node.type || "LOCAL",
        description: "Artifact declared by anvi'o's ProgramsNetwork.",
        requires: [],
        provides: [],
        parameters: [],
        providedByAnvio: Boolean(node.provided_by_anvio)
      });
    }
  });
  const nodes = new Map();
  const edges = [];
  const producersByArtifact = new Map();
  const consumersByArtifact = new Map();
  const essentialInputsByProgram = new Map();
  const programSearchByKey = new Map();

  artifacts.forEach((artifact, id) => nodes.set(`artifact:${id}`, { kind: "artifact", id, item: artifact }));
  programs.forEach((program, key) => {
    essentialInputsByProgram.set(key, baseGraphInputArtifacts(program));
    nodes.set(`program:${key}`, { kind: "program", id: key, item: program });
  });

  const networkNodes = network.nodes || [];
  (network.links || []).forEach((link) => {
    const source = networkNodes[Number(link.source)] || link.source;
    const target = networkNodes[Number(link.target)] || link.target;
    const sourceId = typeof source === "string" ? source : String(source?.id || source?.name || "");
    const targetId = typeof target === "string" ? target : String(target?.id || target?.name || "");
    const sourceIsProgram = (typeof source !== "string" && source?.type === "PROGRAM") || programs.has(sourceId);
    const targetIsProgram = (typeof target !== "string" && target?.type === "PROGRAM") || programs.has(targetId);
    if (sourceIsProgram && !targetIsProgram) {
      const program = programs.get(sourceId);
      const artifact = artifacts.get(targetId.toLowerCase());
      if (!program || !artifact) return;
      edges.push({ from: `program:${sourceId}`, to: `artifact:${artifact.id}`, artifact: artifact.id, kind: "output" });
      if (!producersByArtifact.has(artifact.id)) producersByArtifact.set(artifact.id, []);
      producersByArtifact.get(artifact.id).push(program);
      return;
    }
    if (!sourceIsProgram && targetIsProgram) {
      const artifact = artifacts.get(sourceId.toLowerCase());
      const program = programs.get(targetId);
      if (!artifact || !program) return;
      edges.push({ from: `artifact:${artifact.id}`, to: `program:${targetId}`, artifact: artifact.id, kind: "input" });
      if (!consumersByArtifact.has(artifact.id)) consumersByArtifact.set(artifact.id, []);
      consumersByArtifact.get(artifact.id).push(program);
    }
  });

  return finalizeIntentGraph({
    nodes,
    edges,
    artifacts,
    programs,
    producersByArtifact,
    consumersByArtifact,
    essentialInputsByProgram,
    programSearchByKey,
    source: "anvio-programs-network"
  });
}

function buildIntentGraphFromRuntimeMetadata(version) {
  const artifacts = new Map((version.dataBlocks || []).map((artifact) => [artifact.id, artifact]));
  const programs = new Map(
    (version.programs || [])
      .filter(isInternalProgramItem)
      .map((program) => [programItemKey(program), program])
  );
  programs.forEach((program) => {
    [...(program.requires || []), ...(program.provides || [])].forEach((artifactId) => {
      if (artifactId && !artifacts.has(artifactId)) artifacts.set(artifactId, findDataBlock(artifactId));
    });
  });
  const nodes = new Map();
  const edges = [];
  const producersByArtifact = new Map();
  const consumersByArtifact = new Map();
  const essentialInputsByProgram = new Map();
  const programSearchByKey = new Map();

  artifacts.forEach((artifact, id) => nodes.set(`artifact:${id}`, { kind: "artifact", id, item: artifact }));
  programs.forEach((program, key) => {
    essentialInputsByProgram.set(key, baseGraphInputArtifacts(program));
    nodes.set(`program:${key}`, { kind: "program", id: key, item: program });
    (program.requires || []).forEach((artifactId) => {
      const artifact = artifacts.get(artifactId) || findDataBlock(artifactId);
      edges.push({ from: `artifact:${artifact.id}`, to: `program:${key}`, artifact: artifact.id, kind: "input" });
      if (!consumersByArtifact.has(artifact.id)) consumersByArtifact.set(artifact.id, []);
      consumersByArtifact.get(artifact.id).push(program);
    });
    (program.provides || []).forEach((artifactId) => {
      const artifact = artifacts.get(artifactId) || findDataBlock(artifactId);
      edges.push({ from: `program:${key}`, to: `artifact:${artifact.id}`, artifact: artifact.id, kind: "output" });
      if (!producersByArtifact.has(artifact.id)) producersByArtifact.set(artifact.id, []);
      producersByArtifact.get(artifact.id).push(program);
    });
  });

  return finalizeIntentGraph({
    nodes,
    edges,
    artifacts,
    programs,
    producersByArtifact,
    consumersByArtifact,
    essentialInputsByProgram,
    programSearchByKey,
    source: "runtime-metadata-fallback"
  });
}

function finalizeIntentGraph(graph) {
  const {
    nodes,
    edges,
    artifacts,
    programs,
    producersByArtifact,
    consumersByArtifact,
    essentialInputsByProgram,
    programSearchByKey,
    source
  } = graph;

  const incomingByNode = new Map([...nodes.keys()].map((id) => [id, []]));
  const outgoingByNode = new Map([...nodes.keys()].map((id) => [id, []]));
  edges.forEach((edge) => {
    if (!outgoingByNode.has(edge.from)) outgoingByNode.set(edge.from, []);
    if (!incomingByNode.has(edge.to)) incomingByNode.set(edge.to, []);
    outgoingByNode.get(edge.from).push(edge);
    incomingByNode.get(edge.to).push(edge);
  });

  const programSearch = [...programs.values()].map((program) => {
    const keywordIndex = graphWeightedKeywordIndex([
      { text: commandKeywordText(program.command || program.id), weight: 18 },
      { text: program.title, weight: 8 },
      { text: program.description, weight: 3.5 },
      { text: program.documentation, weight: 1 },
      { text: (program.requires || []).join(" "), weight: 11 },
      { text: (program.provides || []).join(" "), weight: 12 },
      { text: (program.parameters || []).map((parameter) => `${parameter.name || ""} ${(parameter.aliases || []).join(" ")}`).join(" "), weight: 5 },
      { text: (program.parameters || []).map((parameter) => parameter.help || "").join(" "), weight: 1.8 }
    ]);
    const vector = keywordIndex.vector;
    const doc = {
      program,
      key: programItemKey(program),
      nodeId: `program:${programItemKey(program)}`,
      keywords: keywordIndex.keywords,
      keywordSet: keywordIndex.keywordSet,
      keywordWeights: keywordIndex.keywordWeights,
      searchText: keywordIndex.searchText,
      vector,
      magnitude: vectorMagnitude(vector)
    };
    programSearchByKey.set(doc.key, doc);
    return doc;
  });

  const artifactSearch = [...artifacts.values()].map((artifact) => {
    const producers = producersByArtifact.get(artifact.id) || [];
    const consumers = consumersByArtifact.get(artifact.id) || [];
    const keywordIndex = graphWeightedKeywordIndex([
      { text: artifactKeywordText(artifact), weight: 20 },
      { text: artifact.title, weight: 10 },
      { text: artifact.artifactType, weight: 8 },
      { text: artifact.description, weight: 4 },
      { text: producers.map((program) => commandKeywordText(program.command || program.id)).join(" "), weight: 7 },
      { text: consumers.map((program) => commandKeywordText(program.command || program.id)).join(" "), weight: 5 },
      { text: [...(artifact.providedBy || []), ...(artifact.requiredBy || [])].join(" "), weight: 5 }
    ]);
    const vector = keywordIndex.vector;
    return {
      artifact,
      nodeId: `artifact:${artifact.id}`,
      keywords: keywordIndex.keywords,
      keywordSet: keywordIndex.keywordSet,
      keywordWeights: keywordIndex.keywordWeights,
      searchText: keywordIndex.searchText,
      vector,
      magnitude: vectorMagnitude(vector),
      mentionKeys: [
        artifact.id,
        artifact.title,
        ...(artifact.aliases || [])
      ].map(normalizedIntentText).filter((candidate) => candidate.length >= 4)
    };
  });

  return {
    nodes,
    edges,
    artifacts,
    programs,
    producersByArtifact,
    consumersByArtifact,
    essentialInputsByProgram,
    incomingByNode,
    outgoingByNode,
    programSearch,
    programSearchByKey,
    artifactSearch,
    source
  };
}

function graphWeightedKeywordIndex(parts) {
  const keywordWeights = new Map();
  const searchParts = [];
  parts.forEach((part) => {
    const text = typeof part === "string" ? part : part?.text;
    const weight = typeof part === "string" ? 1 : Number(part?.weight || 1);
    if (!text) return;
    searchParts.push(text);
    const tokens = new Set(tokenizeConceptText(text));
    tokens.forEach((token) => {
      keywordWeights.set(token, (keywordWeights.get(token) || 0) + weight);
    });
  });
  const keywords = [...keywordWeights.entries()]
    .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
    .slice(0, 80)
    .map(([token]) => token);
  return {
    keywords,
    keywordSet: new Set(keywords),
    keywordWeights,
    vector: new Map(keywordWeights),
    topKeywords: keywords.slice(0, 18),
    searchText: normalizedIntentText(searchParts.join(" "))
  };
}

function commandKeywordText(command) {
  const value = String(command || "");
  const stripped = value.replace(/^anvi-/, "");
  return [
    value,
    stripped,
    value.replace(/-/g, " "),
    stripped.replace(/-/g, " ")
  ].join(" ");
}

function artifactKeywordText(artifact) {
  const id = String(artifact?.id || "");
  const dbExpanded = id
    .replace(/-db$/i, " database db")
    .replace(/-file$/i, " file")
    .replace(/-files$/i, " files")
    .replace(/-/g, " ");
  return [
    id,
    dbExpanded,
    artifact?.title || "",
    artifact?.artifactType || "",
    ...(artifact?.aliases || [])
  ].join(" ");
}

function graphKeywordIndex(text) {
  const weighted = graphWeightedKeywordIndex([{ text, weight: 1 }]);
  return {
    keywords: weighted.keywords,
    keywordSet: weighted.keywordSet,
    keywordWeights: weighted.keywordWeights,
    vector: weighted.vector,
    searchText: normalizedIntentText(text)
  };
}

function isInternalProgramItem(item) {
  const command = String(item?.command || item?.id || "");
  return item?.kind === "program" && command.startsWith("anvi-") && !item.externalFromAnvio;
}

function artifactCompatible(expected, provided) {
  const left = normalizedDataKey(expected);
  const right = normalizedDataKey(provided);
  if (!left || !right) return false;
  if (left === right) return true;
  if (right.endsWith(`-${left}`)) return true;
  if (left.endsWith(`-${right}`)) return true;
  if (left === "profile-db" && /(^|-)profile-db$/.test(right)) return true;
  if (left === "contigs-db" && /(^|-)contigs-db$/.test(right)) return true;
  if (left === "bam-file" && /(^|-)bam-file$/.test(right)) return true;
  return false;
}

function hasCompatibleArtifactInSet(artifactId, artifactSet) {
  return [...artifactSet].some((candidate) => artifactCompatible(artifactId, candidate));
}

function estimateProgramGraphCost(program, queryOrContext, plannedPrograms = new Map(), graph = cachedIntentGraph(), seen = new Set()) {
  const context = typeof queryOrContext === "string" ? intentQueryContext(queryOrContext, graph) : queryOrContext;
  const key = programItemKey(program);
  if (!key || seen.has(key)) return 8;
  if (plannedPrograms.has(key) || activeNodes().some((node) => node.kind === "program" && programNodeKey(node) === key)) return 0.4;
  seen.add(key);
  const cost = 1 + graphInputArtifacts(program, context, graph)
    .reduce((sum, artifactId) => sum + estimateArtifactGraphCost(artifactId, context, plannedPrograms, graph, seen), 0);
  seen.delete(key);
  return cost;
}

function estimateArtifactGraphCost(artifactId, queryOrContext, plannedPrograms, graph, seen) {
  const context = typeof queryOrContext === "string" ? intentQueryContext(queryOrContext, graph) : queryOrContext;
  if (existingProducerForArtifact(artifactId)) return 0.2;
  if (activeNodes().some((node) => node.kind === "data" && artifactCompatible(artifactId, node.refId))) return 0.2;
  const producers = graph.producersByArtifact.get(artifactId) || [];
  if (!producers.length) return 1.6;
  return Math.min(...producers.map((producer) => 0.7 + estimateProgramGraphCost(producer, context, plannedPrograms, graph, seen)));
}

function intentGoalTitle(goal) {
  if (goal?.type === "artifact") return goal.artifact.title || goal.artifact.id;
  return goal?.program?.command || goal?.program?.title || "Workflow";
}

function activeGraphCenter(fallback = viewportWorldCenter()) {
  const nodes = activeNodes();
  if (!nodes.length) return fallback;
  const minX = Math.min(...nodes.map((node) => node.x));
  const minY = Math.min(...nodes.map((node) => node.y));
  const maxX = Math.max(...nodes.map((node) => node.x + node.w));
  const maxY = Math.max(...nodes.map((node) => node.y + node.h));
  return {
    x: (minX + maxX) / 2,
    y: (minY + maxY) / 2
  };
}

function applyOptimizedPlacement(center = viewportWorldCenter()) {
  const nodes = activeNodes();
  if (!nodes.length) return;
  normalizeActiveNodeSizes();
  const result = optimizePlacementWithPfgo(nodes, activeEdges(), center);
  if (!result) return;
  applyPfgoPlacementResult(result, center);
  state.selectedEdgeId = "";
}

function optimizePlacementWithPfgo(nodes, tabEdges, center) {
  const nodeIds = nodes.map((node) => node.id);
  const nodeIdSet = new Set(nodeIds);
  const edges = tabEdges
    .filter((edge) => nodeIdSet.has(edge.from) && nodeIdSet.has(edge.to) && edge.from !== edge.to)
    .map((edge) => [edge.from, edge.to, edge.id]);
  if (!nodeIds.length) return null;

  const grid = pfgoGridDimensions(nodes, edges);
  const cells = pfgoCells(grid.width, grid.height);
  const layoutLimit = nodeIds.length <= 10 ? 900 : nodeIds.length <= 18 ? 520 : 260;
  const layouts = pfgoGenerateLayouts(nodeIds, edges, cells, grid.width, grid.height, layoutLimit);
  let best = null;

  layouts.forEach((positions) => {
    const routed = pfgoRouteEdges(positions, edges, grid.width, grid.height);
    const totalLength = routed ? routed.totalLength : pfgoLayoutScore(positions, edges)[0] + edges.length * 4;
    const candidate = {
      positions,
      paths: routed?.paths || new Map(),
      score: pfgoResultScore(positions, edges, totalLength)
    };
    if (!best || pfgoCompareScore(candidate.score, best.score) < 0) best = candidate;
  });

  if (!best) {
    const positions = pfgoFallbackLayeredPositions(nodes, edges, grid.width, grid.height);
    const routed = pfgoRouteEdges(positions, edges, grid.width, grid.height);
    best = {
      positions,
      paths: routed?.paths || new Map(),
      score: pfgoResultScore(positions, edges, routed?.totalLength || pfgoLayoutScore(positions, edges)[0])
    };
  }

  return {
    ...best,
    nodes,
    edges,
    grid,
    center
  };
}

function pfgoGridDimensions(nodes, edges) {
  const layers = pfgoTopologicalLayers(nodes.map((node) => node.id), edges);
  const layerCount = Math.max(1, new Set([...layers.values()]).size);
  let width = clamp(Math.max(layerCount, Math.ceil(Math.sqrt(nodes.length * 1.45))), 3, Math.max(3, nodes.length + 2));
  let height = clamp(Math.ceil(nodes.length / width) + 2, 3, Math.max(3, nodes.length + 2));
  while (width * height < nodes.length) width += 1;
  return { width, height };
}

function pfgoCells(width, height) {
  const cells = [];
  for (let y = 0; y < height; y += 1) {
    for (let x = 0; x < width; x += 1) cells.push({ x, y });
  }
  return cells;
}

function pfgoGenerateLayouts(nodes, edges, cells, width, height, limit) {
  const seen = new Set();
  const layouts = [];
  const addLayout = (layout) => {
    const key = nodes.map((node) => `${node}:${pfgoCellKey(layout.get(node))}`).sort().join("|");
    if (seen.has(key)) return;
    seen.add(key);
    layouts.push(layout);
  };

  pfgoGenerateBlockLayouts(nodes, edges, cells, width, height, Math.max(80, Math.floor(limit * 0.7))).forEach(addLayout);
  pfgoGenerateRankedLayouts(nodes, edges, cells, width, height, Math.max(40, limit - layouts.length)).forEach(addLayout);
  if (!layouts.length) addLayout(pfgoFallbackLayeredPositions(nodes.map((id) => ({ id })), edges, width, height));
  return layouts.slice(0, limit);
}

function pfgoGenerateBlockLayouts(nodes, edges, cells, width, height, limit) {
  if (!edges.length || limit <= 0) return [];
  const degrees = pfgoNodeDegrees(nodes, edges);
  const adjacency = pfgoUndirectedAdjacency(nodes, edges);
  const seedEdge = [...edges].sort((a, b) => (degrees.get(b[0]) + degrees.get(b[1])) - (degrees.get(a[0]) + degrees.get(a[1])) || a[0].localeCompare(b[0]) || a[1].localeCompare(b[1]))[0];
  const startNodes = [seedEdge[0], seedEdge[1]];
  let states = [];

  pfgoAdjacentCellPairs(cells, width, height).slice(0, Math.min(48, cells.length * 2)).forEach(([first, second]) => {
    const placed = [first, second];
    states.push({
      placedNodes: startNodes,
      placedCells: placed,
      remainingCells: cells.filter((cell) => !pfgoSameCell(cell, first) && !pfgoSameCell(cell, second)),
      score: pfgoManhattan(first, second),
      lower: pfgoLayoutLowerBound(startNodes, placed, nodes.filter((node) => !startNodes.includes(node)), cells, edges)
    });
  });
  states = pfgoSelectBlockStates(states, limit);

  while (states.length && states[0].placedNodes.length < nodes.length) {
    const expanded = [];
    states.forEach((stateItem) => {
      const assigned = pfgoAssignedMap(stateItem.placedNodes, stateItem.placedCells);
      const remainingNodes = nodes.filter((node) => !assigned.has(node));
      const nextNodes = pfgoNextBlockNodes(remainingNodes, assigned, adjacency, degrees).slice(0, 3);
      nextNodes.forEach((node) => {
        pfgoCandidateBlockCells(node, assigned, stateItem.remainingCells, adjacency, width, height).slice(0, 12).forEach((cell) => {
          const placedNodes = [...stateItem.placedNodes, node];
          const placedCells = [...stateItem.placedCells, cell];
          const remainingCells = stateItem.remainingCells.filter((candidate) => !pfgoSameCell(candidate, cell));
          const remainingAfter = nodes.filter((candidate) => !placedNodes.includes(candidate));
          expanded.push({
            placedNodes,
            placedCells,
            remainingCells,
            score: stateItem.score + pfgoNewEdgeScore(node, cell, assigned, edges),
            lower: pfgoLayoutLowerBound(placedNodes, placedCells, remainingAfter, remainingCells, edges)
          });
        });
      });
    });
    states = pfgoSelectBlockStates(expanded, limit);
  }

  return states
    .filter((stateItem) => stateItem.placedNodes.length === nodes.length)
    .map((stateItem) => pfgoAssignedMap(stateItem.placedNodes, stateItem.placedCells))
    .sort((a, b) => pfgoCompareScore(pfgoLayoutScore(a, edges), pfgoLayoutScore(b, edges)))
    .slice(0, limit);
}

function pfgoGenerateRankedLayouts(nodes, edges, cells, width, height, limit) {
  const orderedNodes = pfgoOrderNodesByConnectivity(nodes, edges);
  const orderedCells = [...cells].sort((a, b) => pfgoCenterDistance(a, width, height) - pfgoCenterDistance(b, width, height) || pfgoCellKey(a).localeCompare(pfgoCellKey(b)));
  let states = [{ placed: [], remainingCells: orderedCells, score: 0, centerScore: 0, lower: 0 }];

  orderedNodes.forEach((node, index) => {
    const assignedNodes = orderedNodes.slice(0, index);
    const remainingNodes = orderedNodes.slice(index + 1);
    const expanded = [];
    states.forEach((stateItem) => {
      const assigned = pfgoAssignedMap(assignedNodes, stateItem.placed);
      stateItem.remainingCells.forEach((cell, cellIndex) => {
        const placed = [...stateItem.placed, cell];
        const remainingCells = stateItem.remainingCells.filter((_, indexCandidate) => indexCandidate !== cellIndex);
        expanded.push({
          placed,
          remainingCells,
          score: stateItem.score + pfgoNewEdgeScore(node, cell, assigned, edges),
          centerScore: stateItem.centerScore + pfgoCenterDistance(cell, width, height),
          lower: pfgoLayoutLowerBound([...assignedNodes, node], placed, remainingNodes, remainingCells, edges)
        });
      });
    });
    states = pfgoSelectLayoutStates(expanded, limit);
  });

  return states.map((stateItem) => pfgoAssignedMap(orderedNodes, stateItem.placed));
}

function pfgoSelectLayoutStates(states, limit) {
  return states
    .sort((a, b) => a.lower - b.lower || a.score - b.score || a.centerScore - b.centerScore || pfgoCellListKey(a.placed).localeCompare(pfgoCellListKey(b.placed)))
    .slice(0, limit);
}

function pfgoSelectBlockStates(states, limit) {
  const compact = states
    .sort((a, b) => a.lower - b.lower || a.score - b.score || pfgoCellListKey(a.placedCells).localeCompare(pfgoCellListKey(b.placedCells)))
    .slice(0, Math.max(1, Math.floor(limit * 0.65)));
  const diverse = states
    .sort((a, b) => pfgoSpreadScore(b.placedCells) - pfgoSpreadScore(a.placedCells) || a.score - b.score)
    .slice(0, Math.max(1, limit - compact.length));
  const seen = new Set();
  return [...compact, ...diverse].filter((stateItem) => {
    const key = `${stateItem.placedNodes.join(",")}:${pfgoCellListKey(stateItem.placedCells)}`;
    if (seen.has(key)) return false;
    seen.add(key);
    return true;
  }).slice(0, limit);
}

function pfgoRouteEdges(positions, edges, width, height) {
  const nodeCells = new Map([...positions.entries()].map(([node, cell]) => [pfgoCellKey(cell), node]));
  const degrees = pfgoNodeDegrees([...positions.keys()], edges);
  const orderedEdges = [...edges].sort((a, b) => pfgoManhattan(positions.get(a[0]), positions.get(a[1])) - pfgoManhattan(positions.get(b[0]), positions.get(b[1]))
    || (degrees.get(b[0]) + degrees.get(b[1])) - (degrees.get(a[0]) + degrees.get(a[1])));
  const routed = new Map();
  let totalLength = 0;

  for (const edge of orderedEdges) {
    const path = pfgoFindShortestPath(positions.get(edge[0]), positions.get(edge[1]), nodeCells, routed, edge, width, height);
    if (!path) return null;
    routed.set(edge[2], path);
    totalLength += path.length - 1;
  }
  return { paths: routed, totalLength };
}

function pfgoFindShortestPath(source, target, nodeCells, routedPaths, currentEdge, width, height) {
  const queue = [source];
  const distance = new Map([[pfgoCellKey(source), 0]]);
  const parents = new Map([[pfgoCellKey(source), []]]);
  let bestDistance = null;

  while (queue.length) {
    const position = queue.shift();
    const step = distance.get(pfgoCellKey(position));
    if (bestDistance != null && step >= bestDistance) continue;
    pfgoNeighbors(position, width, height).forEach((neighbor) => {
      if (!pfgoMoveAllowed(neighbor, target, nodeCells, routedPaths, currentEdge)) return;
      const key = pfgoCellKey(neighbor);
      const newDistance = step + 1;
      const known = distance.get(key);
      if (known == null) {
        distance.set(key, newDistance);
        parents.set(key, [position]);
        if (pfgoSameCell(neighbor, target)) bestDistance = newDistance;
        else queue.push(neighbor);
      } else if (known === newDistance) {
        parents.get(key).push(position);
      }
    });
  }

  if (!distance.has(pfgoCellKey(target))) return null;
  const paths = [];
  const build = (cursor, suffix) => {
    if (paths.length >= 16) return;
    if (pfgoSameCell(cursor, source)) {
      paths.push([source, ...suffix]);
      return;
    }
    [...(parents.get(pfgoCellKey(cursor)) || [])]
      .sort((a, b) => pfgoCellKey(a).localeCompare(pfgoCellKey(b)))
      .forEach((parent) => build(parent, [cursor, ...suffix]));
  };
  build(target, []);
  return paths.sort((a, b) => pfgoPathScore(a)[0] - pfgoPathScore(b)[0] || pfgoPathScore(a)[1] - pfgoPathScore(b)[1])[0] || null;
}

function pfgoMoveAllowed(position, target, nodeCells, routedPaths, currentEdge) {
  const key = pfgoCellKey(position);
  if (nodeCells.has(key) && !pfgoSameCell(position, target)) return false;
  for (const [edgeId, path] of routedPaths.entries()) {
    if (!path.some((cell) => pfgoSameCell(cell, position))) continue;
    if (!pfgoSharedEndpointCell(currentEdge, edgeId, position, nodeCells)) return false;
  }
  return true;
}

function pfgoSharedEndpointCell(edge, routedEdgeId, position, nodeCells) {
  const node = nodeCells.get(pfgoCellKey(position));
  if (!node) return false;
  const other = activeEdges().find((candidate) => candidate.id === routedEdgeId);
  if (!other) return false;
  return [edge[0], edge[1]].includes(node) && [other.from, other.to].includes(node);
}

function applyPfgoPlacementResult(result, center) {
  const nodes = result.nodes;
  const edges = activeEdges();
  const maxWidth = Math.max(...nodes.map((node) => node.w || NODE_DEFAULT_WIDTH), NODE_DEFAULT_WIDTH);
  const maxHeight = Math.max(...nodes.map((node) => node.h || NODE_DEFAULT_HEIGHT), NODE_DEFAULT_HEIGHT);
  const cellWidth = Math.max(300, maxWidth + 110);
  const cellHeight = Math.max(210, maxHeight + 78);
  const totalWidth = result.grid.width * cellWidth;
  const totalHeight = result.grid.height * cellHeight;
  const startX = center.x - totalWidth / 2 + cellWidth / 2;
  const startY = center.y - totalHeight / 2 + cellHeight / 2;
  const nodeByLayoutId = new Map(nodes.map((node) => [node.id, node]));
  const cellToPoint = (cell) => ({
    x: Math.round(startX + cell.x * cellWidth),
    y: Math.round(startY + cell.y * cellHeight)
  });

  nodes.forEach((node) => {
    const cell = result.positions.get(node.id);
    if (!cell) return;
    const point = cellToPoint(cell);
    node.x = clamp(Math.round(point.x - node.w / 2), 20, WORLD_WIDTH - node.w - 20);
    node.y = clamp(Math.round(point.y - node.h / 2), 20, WORLD_HEIGHT - node.h - 20);
  });

  edges.forEach((edge) => {
    edge.points = [];
    delete edge.layoutRoute;
    const path = result.paths.get(edge.id);
    const from = nodeByLayoutId.get(edge.from);
    const to = nodeByLayoutId.get(edge.to);
    if (!path || !from || !to || path.length < 2) return;
    edge.layoutRoute = pfgoPathToCanvasRoute(path, from, to, cellToPoint);
  });
}

function pfgoPathToCanvasRoute(path, from, to, cellToPoint) {
  const firstDirection = pfgoDirection(path[0], path[1]);
  const lastDirection = pfgoDirection(path[path.length - 2], path[path.length - 1]);
  const start = pfgoNodePort(from, firstDirection);
  const end = pfgoNodePort(to, { x: -lastDirection.x, y: -lastDirection.y });
  const route = [start];
  path.slice(1, -1).map(cellToPoint).forEach((point) => pfgoAppendOrthogonalPoint(route, point));
  pfgoAppendOrthogonalPoint(route, end);
  return simplifyRoutePoints(route);
}

function pfgoNodePort(node, direction) {
  const cx = node.x + node.w / 2;
  const cy = node.y + node.h / 2;
  if (Math.abs(direction.x) >= Math.abs(direction.y) && direction.x > 0) return { x: node.x + node.w, y: cy };
  if (Math.abs(direction.x) >= Math.abs(direction.y) && direction.x < 0) return { x: node.x, y: cy };
  if (direction.y > 0) return { x: cx, y: node.y + node.h };
  if (direction.y < 0) return { x: cx, y: node.y };
  return { x: node.x + node.w, y: cy };
}

function pfgoAppendOrthogonalPoint(route, point) {
  const current = route[route.length - 1];
  if (!current) {
    route.push(point);
    return;
  }
  if (current.x !== point.x && current.y !== point.y) route.push({ x: point.x, y: current.y });
  route.push(point);
}

function simplifyRoutePoints(points) {
  const rounded = points.map((point) => ({ x: Math.round(point.x), y: Math.round(point.y) }));
  const compact = [];
  rounded.forEach((point) => {
    const previous = compact[compact.length - 1];
    if (!previous || previous.x !== point.x || previous.y !== point.y) compact.push(point);
  });
  return compact.filter((point, index, array) => {
    if (index === 0 || index === array.length - 1) return true;
    const previous = array[index - 1];
    const next = array[index + 1];
    return !((previous.x === point.x && point.x === next.x) || (previous.y === point.y && point.y === next.y));
  });
}

function pfgoFallbackLayeredPositions(nodes, edges, width, height) {
  const ids = nodes.map((node) => typeof node === "string" ? node : node.id);
  const layers = pfgoTopologicalLayers(ids, edges);
  const groups = new Map();
  ids.forEach((id) => {
    const layer = layers.get(id) || 0;
    if (!groups.has(layer)) groups.set(layer, []);
    groups.get(layer).push(id);
  });
  const positions = new Map();
  [...groups.entries()].sort((a, b) => a[0] - b[0]).forEach(([layer, group]) => {
    group.sort().forEach((id, index) => {
      positions.set(id, {
        x: clamp(layer, 0, width - 1),
        y: clamp(index + Math.floor((height - group.length) / 2), 0, height - 1)
      });
    });
  });
  return positions;
}

function pfgoTopologicalLayers(nodes, edges) {
  const ids = nodes.map((node) => typeof node === "string" ? node : node.id);
  const outgoing = new Map(ids.map((id) => [id, []]));
  const indegree = new Map(ids.map((id) => [id, 0]));
  edges.forEach(([from, to]) => {
    if (!outgoing.has(from) || !outgoing.has(to)) return;
    outgoing.get(from).push(to);
    indegree.set(to, (indegree.get(to) || 0) + 1);
  });
  const layers = new Map(ids.map((id) => [id, 0]));
  const queue = ids.filter((id) => (indegree.get(id) || 0) === 0).sort();
  while (queue.length) {
    const id = queue.shift();
    outgoing.get(id).forEach((to) => {
      layers.set(to, Math.max(layers.get(to) || 0, (layers.get(id) || 0) + 1));
      indegree.set(to, (indegree.get(to) || 0) - 1);
      if ((indegree.get(to) || 0) === 0) queue.push(to);
    });
    queue.sort();
  }
  return layers;
}

function pfgoAssignedMap(nodes, cells) {
  return new Map(nodes.map((node, index) => [node, cells[index]]));
}

function pfgoNodeDegrees(nodes, edges) {
  const degrees = new Map(nodes.map((node) => [node, 0]));
  edges.forEach(([from, to]) => {
    degrees.set(from, (degrees.get(from) || 0) + 1);
    degrees.set(to, (degrees.get(to) || 0) + 1);
  });
  return degrees;
}

function pfgoUndirectedAdjacency(nodes, edges) {
  const adjacency = new Map(nodes.map((node) => [node, new Set()]));
  edges.forEach(([from, to]) => {
    adjacency.get(from)?.add(to);
    adjacency.get(to)?.add(from);
  });
  return adjacency;
}

function pfgoNextBlockNodes(remainingNodes, assigned, adjacency, degrees) {
  const scored = remainingNodes.map((node) => ({
    node,
    placedNeighbors: [...(adjacency.get(node) || [])].filter((neighbor) => assigned.has(neighbor)).length,
    degree: degrees.get(node) || 0
  }));
  const frontier = Math.max(...scored.map((item) => item.placedNeighbors), 0);
  return scored
    .filter((item) => frontier <= 0 || item.placedNeighbors === frontier)
    .sort((a, b) => b.placedNeighbors - a.placedNeighbors || b.degree - a.degree || a.node.localeCompare(b.node))
    .map((item) => item.node);
}

function pfgoCandidateBlockCells(node, assigned, remainingCells, adjacency, width, height) {
  const neighborPositions = [...(adjacency.get(node) || [])].filter((neighbor) => assigned.has(neighbor)).map((neighbor) => assigned.get(neighbor));
  return [...remainingCells].sort((a, b) => {
    const aScore = neighborPositions.length ? neighborPositions.reduce((sum, cell) => sum + pfgoManhattan(a, cell), 0) : pfgoCenterDistance(a, width, height);
    const bScore = neighborPositions.length ? neighborPositions.reduce((sum, cell) => sum + pfgoManhattan(b, cell), 0) : pfgoCenterDistance(b, width, height);
    return aScore - bScore
      || pfgoCenterDistance(a, width, height) - pfgoCenterDistance(b, width, height)
      || pfgoCellKey(a).localeCompare(pfgoCellKey(b));
  });
}

function pfgoOrderNodesByConnectivity(nodes, edges) {
  const degrees = pfgoNodeDegrees(nodes, edges);
  return [...nodes].sort((a, b) => (degrees.get(b) || 0) - (degrees.get(a) || 0) || a.localeCompare(b));
}

function pfgoAdjacentCellPairs(cells, width, height) {
  const pairs = [];
  cells.forEach((cell) => {
    pfgoNeighbors(cell, width, height).forEach((neighbor) => pairs.push([cell, neighbor]));
  });
  return pairs.sort((a, b) => pfgoCenterDistance(a[0], width, height) + pfgoCenterDistance(a[1], width, height)
    - pfgoCenterDistance(b[0], width, height) - pfgoCenterDistance(b[1], width, height));
}

function pfgoNeighbors(cell, width, height) {
  return [
    { x: cell.x + 1, y: cell.y },
    { x: cell.x - 1, y: cell.y },
    { x: cell.x, y: cell.y + 1 },
    { x: cell.x, y: cell.y - 1 }
  ].filter((candidate) => candidate.x >= 0 && candidate.x < width && candidate.y >= 0 && candidate.y < height);
}

function pfgoNewEdgeScore(node, cell, assigned, edges) {
  return edges.reduce((score, [from, to]) => {
    if (from === node && assigned.has(to)) return score + pfgoManhattan(cell, assigned.get(to));
    if (to === node && assigned.has(from)) return score + pfgoManhattan(cell, assigned.get(from));
    return score;
  }, 0);
}

function pfgoLayoutLowerBound(assignedNodes, assignedCells, remainingNodes, remainingCells, edges) {
  const assigned = pfgoAssignedMap(assignedNodes, assignedCells);
  const remaining = new Set(remainingNodes);
  return edges.reduce((total, [from, to]) => {
    const fromCell = assigned.get(from);
    const toCell = assigned.get(to);
    if (fromCell && toCell) return total + pfgoManhattan(fromCell, toCell);
    if (fromCell && remaining.has(to) && remainingCells.length) return total + Math.min(...remainingCells.map((cell) => pfgoManhattan(fromCell, cell)));
    if (toCell && remaining.has(from) && remainingCells.length) return total + Math.min(...remainingCells.map((cell) => pfgoManhattan(toCell, cell)));
    if (remaining.has(from) && remaining.has(to)) return total + 1;
    return total;
  }, 0);
}

function pfgoLayoutScore(positions, edges) {
  const incident = new Map([...positions.keys()].map((node) => [node, 0]));
  let total = 0;
  edges.forEach(([from, to]) => {
    const distance = pfgoManhattan(positions.get(from), positions.get(to));
    total += distance;
    incident.set(from, (incident.get(from) || 0) + distance);
    incident.set(to, (incident.get(to) || 0) + distance);
  });
  return [total, Math.max(...incident.values(), 0), pfgoSpreadScore([...positions.values()])];
}

function pfgoResultScore(positions, edges, totalLength) {
  const incident = new Map([...positions.keys()].map((node) => [node, 0]));
  edges.forEach(([from, to]) => {
    const distance = pfgoManhattan(positions.get(from), positions.get(to));
    incident.set(from, (incident.get(from) || 0) + distance);
    incident.set(to, (incident.get(to) || 0) + distance);
  });
  return [totalLength, Math.max(...incident.values(), 0), ...pfgoLayoutScore(positions, edges)];
}

function pfgoCompareScore(a, b) {
  for (let i = 0; i < Math.max(a.length, b.length); i += 1) {
    const diff = (a[i] || 0) - (b[i] || 0);
    if (diff) return diff;
  }
  return 0;
}

function pfgoPathScore(path) {
  let turns = 0;
  for (let index = 2; index < path.length; index += 1) {
    const a = pfgoDirection(path[index - 2], path[index - 1]);
    const b = pfgoDirection(path[index - 1], path[index]);
    if (a.x !== b.x || a.y !== b.y) turns += 1;
  }
  return [path.length, turns];
}

function pfgoDirection(a, b) {
  return { x: Math.sign(b.x - a.x), y: Math.sign(b.y - a.y) };
}

function pfgoManhattan(a, b) {
  if (!a || !b) return 0;
  return Math.abs(a.x - b.x) + Math.abs(a.y - b.y);
}

function pfgoCenterDistance(cell, width, height) {
  return Math.abs((2 * cell.x) - (width - 1)) + Math.abs((2 * cell.y) - (height - 1));
}

function pfgoSpreadScore(cells) {
  let total = 0;
  cells.forEach((cell, index) => {
    for (let otherIndex = index + 1; otherIndex < cells.length; otherIndex += 1) total += pfgoManhattan(cell, cells[otherIndex]);
  });
  return total;
}

function pfgoCellKey(cell) {
  return `${cell?.x ?? 0},${cell?.y ?? 0}`;
}

function pfgoCellListKey(cells) {
  return cells.map(pfgoCellKey).join("|");
}

function pfgoSameCell(a, b) {
  return Boolean(a && b && a.x === b.x && a.y === b.y);
}

function compareLayoutNodes(a, b) {
  const kindScore = { data: 0, program: 1, workflow: 2 };
  return (kindScore[a?.kind] || 0) - (kindScore[b?.kind] || 0)
    || String(a?.title || "").localeCompare(String(b?.title || ""))
    || String(a?.id || "").localeCompare(String(b?.id || ""));
}

function tightenDataNodeLayers(nodes, incoming, outgoing, layers) {
  nodes
    .filter((node) => node.kind === "data")
    .forEach((node) => {
      const parents = incoming.get(node.id) || [];
      const children = outgoing.get(node.id) || [];
      if (!parents.length && children.length) {
        const firstConsumerLayer = Math.min(...children.map((id) => layers.get(id) || 0));
        layers.set(node.id, Math.max(0, firstConsumerLayer - 1));
        return;
      }
      if (parents.length && !children.length) {
        const lastProducerLayer = Math.max(...parents.map((id) => layers.get(id) || 0));
        layers.set(node.id, lastProducerLayer + 1);
      }
    });
}

function applyImportOrderLayers(nodes, incoming, outgoing, layers) {
  const imported = nodes
    .filter((node) => Number.isFinite(node.importOrder))
    .sort((a, b) => a.importOrder - b.importOrder);
  if (!imported.length) return;
  const connectedIds = new Set();
  imported.forEach((node) => {
    if ((incoming.get(node.id) || []).length || (outgoing.get(node.id) || []).length) connectedIds.add(node.id);
  });
  imported.forEach((node) => {
    const orderLayer = Math.max(0, Math.floor(node.importOrder));
    if (!connectedIds.has(node.id)) {
      layers.set(node.id, orderLayer);
      return;
    }
    layers.set(node.id, Math.max(layers.get(node.id) || 0, Math.min(orderLayer, imported.length - 1)));
  });
}

function orderLayoutGroupsByConnections(orderedLayers, incoming, outgoing) {
  orderedLayers.forEach(([, group]) => group.sort(compareLayoutNodes));
  const order = () => new Map(orderedLayers.flatMap(([, group]) => group.map((node, index) => [node.id, index])));
  for (let pass = 0; pass < 4; pass += 1) {
    let orderMap = order();
    for (let index = 1; index < orderedLayers.length; index += 1) {
      const [, group] = orderedLayers[index];
      group.sort((a, b) => compareByAdjacentOrder(a, b, incoming, orderMap));
    }
    orderMap = order();
    for (let index = orderedLayers.length - 2; index >= 0; index -= 1) {
      const [, group] = orderedLayers[index];
      group.sort((a, b) => compareByAdjacentOrder(a, b, outgoing, orderMap));
    }
  }
}

function compareByAdjacentOrder(a, b, adjacency, orderMap) {
  const aScore = adjacentOrderScore(a.id, adjacency, orderMap);
  const bScore = adjacentOrderScore(b.id, adjacency, orderMap);
  return aScore - bScore || compareLayoutNodes(a, b);
}

function adjacentOrderScore(id, adjacency, orderMap) {
  const linked = (adjacency.get(id) || []).filter((linkedId) => orderMap.has(linkedId));
  if (!linked.length) return Number.MAX_SAFE_INTEGER;
  return average(linked.map((linkedId) => orderMap.get(linkedId)));
}

function assignGridRows(orderedLayers, incoming, outgoing, maxRows) {
  const rows = new Map();
  orderedLayers.forEach(([, group]) => {
    group.forEach((node, rowIndex) => {
      rows.set(node.id, centeredGridRow(rowIndex, group.length, maxRows));
    });
  });
  for (let pass = 0; pass < 8; pass += 1) {
    for (let layerIndex = 0; layerIndex < orderedLayers.length; layerIndex += 1) {
      assignLayerRows(orderedLayers[layerIndex][1], incoming, outgoing, rows, maxRows);
    }
    for (let layerIndex = orderedLayers.length - 1; layerIndex >= 0; layerIndex -= 1) {
      assignLayerRows(orderedLayers[layerIndex][1], incoming, outgoing, rows, maxRows);
    }
  }
  return rows;
}

function centeredGridRow(index, count, maxRows) {
  if (count >= maxRows) return index;
  const start = Math.floor((maxRows - count) / 2);
  return clamp(start + index, 0, maxRows - 1);
}

function assignLayerRows(group, incoming, outgoing, rows, maxRows) {
  const desired = group.map((node, index) => {
    const linkedRows = [
      ...(incoming.get(node.id) || []),
      ...(outgoing.get(node.id) || [])
    ].map((id) => rows.get(id)).filter((row) => Number.isFinite(row));
    return {
      node,
      desired: linkedRows.length ? average(linkedRows) : centeredGridRow(index, group.length, maxRows),
      fallback: centeredGridRow(index, group.length, maxRows)
    };
  });
  const assigned = nearestFreeGridRows(desired, maxRows);
  assigned.forEach((row, id) => rows.set(id, row));
}

function nearestFreeGridRows(items, maxRows) {
  const used = new Set();
  const assigned = new Map();
  items
    .sort((a, b) => a.desired - b.desired || compareLayoutNodes(a.node, b.node))
    .forEach((item) => {
      const preferred = clamp(Math.round(item.desired), 0, maxRows - 1);
      let best = preferred;
      for (let distance = 0; distance < maxRows; distance += 1) {
        const candidates = [preferred - distance, preferred + distance]
          .filter((row, index, array) => row >= 0 && row < maxRows && array.indexOf(row) === index);
        const free = candidates.find((row) => !used.has(row));
        if (free != null) {
          best = free;
          break;
        }
      }
      used.add(best);
      assigned.set(item.node.id, best);
    });
  return assigned;
}

function average(values) {
  if (!values.length) return 0;
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}

function render() {
  renderSourceStatus();
  renderTabs();
  renderSearch();
  normalizeActiveNodeSizes();
  renderNodes();
  renderEdges();
  renderInspector();
  renderToolbarState();
}

function renderSourceStatus() {
  if (!dom.sourceStatus) return;
  const version = currentVersion();
  dom.sourceStatus.textContent = version.label;
}

function renderTabs() {
  dom.tabs.innerHTML = state.tabs.map((tab) => `
    <div class="tab ${tab.id === state.activeTabId ? "active" : ""}" role="button" tabindex="0" data-tab-id="${tab.id}" title="Double-click to rename">
      <span class="tab-name">${escapeHtml(tab.name)}</span>
    </div>
  `).join("");
}

function handleTabDoubleClick(event) {
  const tabButton = event.target.closest("[data-tab-id]");
  if (!tabButton || event.target.closest("[data-tab-name-input]")) return;
  event.preventDefault();
  event.stopPropagation();
  beginTabNameEdit(tabButton.dataset.tabId);
}

function beginTabNameEdit(tabId) {
  const tab = state.tabs.find((candidate) => candidate.id === tabId);
  const tabButton = [...dom.tabs.querySelectorAll("[data-tab-id]")]
    .find((candidate) => candidate.dataset.tabId === tabId);
  const name = tabButton?.querySelector(".tab-name");
  if (!tab || !name || name.querySelector("[data-tab-name-input]")) return;

  const input = document.createElement("input");
  input.type = "text";
  input.value = tab.name;
  input.className = "tab-name-input";
  input.dataset.tabNameInput = tabId;
  input.setAttribute("aria-label", "Canvas name");
  name.textContent = "";
  name.appendChild(input);
  input.focus();
  input.select();
}

function commitTabNameEdit(input) {
  const tab = state.tabs.find((candidate) => candidate.id === input.dataset.tabNameInput);
  if (!tab) return;
  const nextName = input.value.trim() || tab.name;
  if (nextName !== tab.name) {
    pushHistory();
    tab.name = nextName;
  }
  render();
}

function cancelTabNameEdit() {
  renderTabs();
}

function handleTabEditKeydown(event) {
  const input = event.target.closest("[data-tab-name-input]");
  if (!input) {
    const tabButton = event.target.closest("[data-tab-id]");
    if (tabButton && (event.key === "Enter" || event.key === " ")) {
      event.preventDefault();
      state.activeTabId = tabButton.dataset.tabId;
      state.selectedNodeId = "";
      state.selectedEdgeId = "";
      render();
    }
    return;
  }
  event.stopPropagation();
  if (event.key === "Enter") {
    event.preventDefault();
    commitTabNameEdit(input);
  }
  if (event.key === "Escape") {
    event.preventDefault();
    cancelTabNameEdit();
  }
}

function handleTabEditFocusout(event) {
  const input = event.target.closest("[data-tab-name-input]");
  if (!input) return;
  commitTabNameEdit(input);
}

function handleLibrarySearchInput() {
  state.searchQuery = String(dom.searchInput.value || "").trim();
  if (state.searchQuery) state.focusedArtifact = null;
  renderSearch();
}

function renderSearch() {
  const version = currentVersion();
  const query = state.searchQuery.toLowerCase().trim();
  if (!query) {
    if (state.focusedArtifact) {
      renderArtifactSearch(state.focusedArtifact);
      return;
    }
    const selectedNode = activeNodes().find((candidate) => candidate.id === state.selectedNodeId);
    if (selectedNode) {
      renderRelatedSearch(selectedNode);
      return;
    }
    dom.libraryMeta.textContent = `${version.label}: type to search documented blocks.`;
    dom.resultList.innerHTML = "";
    return;
  }
  const words = [...new Set(tokenizeConceptText(query))];
  const items = libraryItems()
    .map((item) => ({ item, score: scoreSearchItem(item, query, words) }))
    .filter((result) => result.score > 0 && (state.searchKind === "all" || result.item.kind === state.searchKind))
    .sort((a, b) => b.score - a.score || searchResultTitle(a.item).localeCompare(searchResultTitle(b.item)))
    .slice(0, 80)
    .map((result) => result.item);

  dom.libraryMeta.textContent = `${version.label}: ${items.length} shown from ${version.programs.length + version.dataBlocks.length + version.workflows.length} documented blocks.`;
  renderSearchCards(items, `<div class="result-card"><strong>No match</strong><p>Try a documented term such as contigs, profile, FASTA, taxonomy, or workflow.</p></div>`);
}

function renderArtifactSearch(focus) {
  const related = artifactSearchResults(focus.artifact, focus.role)
    .filter((result) => state.searchKind === "all" || result.item.kind === state.searchKind);
  const roleText = focus.role === "output"
    ? "consumers"
    : focus.role === "input"
      ? "producers"
      : "related blocks";
  dom.libraryMeta.textContent = `${currentVersion().label}: ${related.length} ${roleText} for ${focus.artifact}.`;
  renderSearchCards(
    related,
    `<div class="result-card"><strong>No related block</strong><p>No documented ${roleText} found for ${escapeHtml(focus.artifact)}.</p></div>`
  );
}

function artifactSearchResults(artifact, role = "both") {
  const related = new Map();
  const currentItems = canvasItemKeys();
  const add = (item, relation, rank) => {
    if (!item) return;
    if (state.searchKind !== "all" && item.kind !== state.searchKind) return;
    if (currentItems.has(itemCanvasKey(item))) return;
    const key = `${item.kind}:${item.id}`;
    const connectionScore = graphConnectionScore(item);
    const existing = related.get(key);
    if (existing) {
      if (!existing.relation.includes(relation)) existing.relation = `${existing.relation}, ${relation}`;
      existing.rank = Math.min(existing.rank, rank);
      existing.connectionScore = Math.max(existing.connectionScore, connectionScore);
      return;
    }
    related.set(key, { item, relation, rank, connectionScore });
  };

  const executables = [...currentVersion().programs, ...currentVersion().workflows];
  if (role === "input" || role === "both") {
    executables.forEach((item) => {
      if ((item.provides || []).includes(artifact)) add(item, `creates ${artifact}`, 1);
    });
  }
  if (role === "output" || role === "both") {
    executables.forEach((item) => {
      if ((item.requires || []).includes(artifact)) add(item, `uses ${artifact}`, 2);
    });
  }
  add(findDataBlock(artifact), `data ${artifact}`, role === "input" ? 3 : 1);

  return [...related.values()]
    .sort((a, b) => a.rank - b.rank || b.connectionScore - a.connectionScore || searchResultTitle(a.item).localeCompare(searchResultTitle(b.item)));
}

function renderSearchCards(results, emptyHtml) {
  if (!results.length) {
    dom.resultList.innerHTML = emptyHtml;
    return;
  }
  dom.resultList.innerHTML = results.map((result) => {
    const item = result.item || result;
    const relation = result.relation || "";
    return `
    <button class="result-card" type="button" data-spawn-kind="${item.kind}" data-spawn-id="${item.id}">
      <strong class="result-command">${escapeHtml(searchResultTitle(item))}</strong>
      ${item.title && item.title !== searchResultTitle(item) ? `<span class="result-title">${escapeHtml(item.title)}</span>` : ""}
      <p>${escapeHtml(compactText(item.description, 180))}</p>
      <span class="chips">
        ${relation ? `<span class="chip relation">${escapeHtml(relation)}</span>` : ""}
        <span class="chip ${item.kind}">${escapeHtml(item.kind)}</span>
        ${isExternalAnvioItem(item) ? `<span class="chip external">external tool</span>` : ""}
        ${item.artifactType ? `<span class="chip data">${escapeHtml(item.artifactType)}</span>` : ""}
        ${(item.requires || []).slice(0, 3).map((name) => `<span class="chip">in: ${escapeHtml(name)}</span>`).join("")}
        ${(item.provides || []).slice(0, 3).map((name) => `<span class="chip">out: ${escapeHtml(name)}</span>`).join("")}
      </span>
    </button>
  `;
  }).join("");
}

function renderRelatedSearch(node) {
  const related = relatedSearchResults(node)
    .filter((result) => state.searchKind === "all" || result.item.kind === state.searchKind);
  const title = nodeDisplayTitle(node);
  dom.libraryMeta.textContent = `${currentVersion().label}: ${related.length} related block${related.length === 1 ? "" : "s"} for ${title}.`;
  renderSearchCards(related, `<div class="result-card"><strong>No related block</strong><p>No exact input/output relation is documented for this block in the local database.</p></div>`);
}

function relatedSearchResults(node) {
  const related = new Map();
  const currentItems = canvasItemKeys();
  const add = (item, relation, rank) => {
    if (!item || item.id === node.refId) return;
    if (state.searchKind !== "all" && item.kind !== state.searchKind) return;
    if (currentItems.has(itemCanvasKey(item))) return;
    const key = `${item.kind}:${item.id}`;
    const connectionScore = graphConnectionScore(item, node);
    const existing = related.get(key);
    if (existing) {
      if (!existing.relation.includes(relation)) existing.relation = `${existing.relation}, ${relation}`;
      existing.rank = Math.min(existing.rank, rank);
      existing.connectionScore = Math.max(existing.connectionScore, connectionScore);
      return;
    }
    related.set(key, { item, relation, rank, connectionScore });
  };

  const executableItems = [...currentVersion().programs, ...currentVersion().workflows];
  if (node.kind === "data") {
    executableItems.forEach((item) => {
      if ((item.provides || []).includes(node.refId)) add(item, `produces ${node.refId}`, 1);
      if ((item.requires || []).includes(node.refId)) add(item, `uses ${node.refId}`, 2);
    });
  } else {
    (node.requires || []).forEach((artifactId) => {
      add(findDataBlock(artifactId), `input ${artifactId}`, 1);
      executableItems.forEach((item) => {
        if (item.id !== node.refId && (item.provides || []).includes(artifactId)) {
          add(item, `produces ${artifactId}`, 2);
        }
      });
    });
    (node.provides || []).forEach((artifactId) => {
      add(findDataBlock(artifactId), `output ${artifactId}`, 3);
      executableItems.forEach((item) => {
        if (item.id !== node.refId && (item.requires || []).includes(artifactId)) {
          add(item, `uses ${artifactId}`, 4);
        }
      });
    });
  }

  return [...related.values()]
    .sort((a, b) => a.rank - b.rank || b.connectionScore - a.connectionScore || searchResultTitle(a.item).localeCompare(searchResultTitle(b.item)));
}

function canvasItemKeys() {
  const keys = new Set();
  activeNodes().forEach((node) => {
    keys.add(itemCanvasKey(node));
    if (node.command) keys.add(`command:${node.command}`);
  });
  return keys;
}

function itemCanvasKey(item) {
  if (!item) return "";
  if (item.kind === "data") return `data:${item.refId || item.id}`;
  return `command:${item.command || item.refId || item.id}`;
}

function graphConnectionScore(item, selectedNode = null) {
  if (!item) return 0;
  const itemRequires = new Set(item.requires || []);
  const itemProvides = new Set(item.provides || []);
  let score = 0;
  activeNodes().forEach((node) => {
    if (selectedNode && node.id === selectedNode.id) return;
    const nodeRequires = new Set(node.requires || []);
    const nodeProvides = new Set(node.provides || []);
    if (node.kind === "data") {
      const artifact = node.refId || node.title;
      if (itemRequires.has(artifact)) score += 4;
      if (itemProvides.has(artifact)) score += 3;
      return;
    }
    itemRequires.forEach((artifact) => {
      if (nodeProvides.has(artifact)) score += 4;
      if (nodeRequires.has(artifact)) score += 1;
    });
    itemProvides.forEach((artifact) => {
      if (nodeRequires.has(artifact)) score += 4;
      if (nodeProvides.has(artifact)) score += 1;
    });
  });
  activeEdges().forEach((edge) => {
    if (edge.label && (itemRequires.has(edge.label) || itemProvides.has(edge.label))) score += 1;
  });
  return score;
}

function documentedPipelineItems() {
  const preloaded = Array.isArray(window.ANVIO_PRELOADED_WORKFLOWS)
    ? window.ANVIO_PRELOADED_WORKFLOWS.map((workflow) => ({
      ...workflow,
      kind: "preloaded-workflow",
      id: `preloaded:${workflow.id}`
    }))
    : [];
  if (preloaded.length) return preloaded;

  const version = currentVersion();
  const pipelines = new Map();
  const addPipeline = (source, workflowName = "") => {
    if (!source) return;
    const cleanWorkflow = workflowName || String(source.id || source.title || "")
      .replace(/^workflow:/, "")
      .replace(/-workflow$/, "");
    if (!cleanWorkflow) return;
    const id = `documented-pipeline:${cleanWorkflow}`;
    if (pipelines.has(id)) return;
    const title = cleanWorkflow.endsWith("workflow") ? cleanWorkflow : `${cleanWorkflow} workflow`;
    pipelines.set(id, {
      id,
      kind: "workflow",
      title,
      command: "anvi-run-workflow",
      description: source.description || `${title} extracted from the local anvi'o runtime database.`,
      requires: ["workflow-config"],
      provides: source.provides || [],
      parameters: [
        {
          name: "--workflow",
          value: cleanWorkflow,
          defaultValue: cleanWorkflow,
          required: true,
          enabled: true,
          source: "local-workflow-list",
          documentation: "Workflow name from the local anvi'o runtime database."
        }
      ],
      examples: [`anvi-run-workflow --workflow ${cleanWorkflow} --config-file config.json`],
      docsUrl: source.docsUrl || version.docsUrl,
      documentation: source.documentation || source.description || ""
    });
  };

  (version.workflows || []).forEach((workflow) => {
    const fromParameter = (workflow.parameters || []).find((param) => param.name === "--workflow")?.value;
    const fromCommand = String(workflow.command || "").match(/--workflow\s+([^\s]+)/)?.[1];
    addPipeline(workflow, fromParameter || fromCommand);
  });

  const workflowArtifact = (version.dataBlocks || []).find((item) => item.id === "workflow");
  extractWorkflowNamesFromDocumentation(workflowArtifact?.documentation || workflowArtifact?.description || "")
    .forEach((workflowName) => addPipeline(workflowArtifact, workflowName));

  (version.dataBlocks || [])
    .filter((item) => item.id !== "workflow" && (item.artifactType === "WORKFLOW" || /(^workflow:|-workflow$)/.test(item.id)))
    .forEach((item) => addPipeline(item));

  return [...pipelines.values()].sort((a, b) => a.title.localeCompare(b.title));
}

function extractWorkflowNamesFromDocumentation(text) {
  const names = [];
  const pattern = /-\s+([A-Za-z0-9][A-Za-z0-9 -]*?)\s+workflow\b/gi;
  let match;
  while ((match = pattern.exec(String(text || "")))) {
    const name = canonicalWorkflowName(match[1]);
    if (name && !names.includes(name)) names.push(name);
  }
  return names;
}

function canonicalWorkflowName(name) {
  const compact = String(name || "")
    .trim()
    .replace(/\s+/g, "-")
    .toLowerCase();
  const aliases = {
    "ecophylo": "ecophylo",
    "eco-phylo": "ecophylo",
    "sra-download": "sra-download"
  };
  return aliases[compact] || compact;
}

function findDocumentedPipeline(id) {
  return documentedPipelineItems().find((item) => item.id === id) || null;
}

function workflowConceptItems() {
  return Array.isArray(window.ANVIO_WORKFLOW_CONCEPTS) ? window.ANVIO_WORKFLOW_CONCEPTS : [];
}

function showWorkflowConcepts() {
  const concepts = workflowConceptItems();
  showModal("Workflow concept search", `
    <div class="concept-dialog">
      <label class="field-stack" for="conceptSearchInput">
        <span>Describe what you want to do</span>
        <textarea id="conceptSearchInput" rows="4" placeholder="Example: I have metagenomic reads and I want to assemble, map reads, profile coverage, and summarize bins."></textarea>
      </label>
      <div class="concept-search-actions">
        <button class="secondary-button" type="button" data-action="run-concept-search">Find matching concept</button>
        <span class="concept-count">${concepts.length} concept document${concepts.length === 1 ? "" : "s"} indexed</span>
      </div>
      <div id="conceptSearchResults" class="concept-results"></div>
      <article id="conceptDocument" class="concept-document">
        <p class="concept-placeholder">Describe an analysis goal to find the closest workflow concept.</p>
      </article>
    </div>
  `);
}

function handleModalInput(event) {
  if (event.target?.id !== "conceptSearchInput") return;
  runConceptSearchFromDialog();
}

function runConceptSearchFromDialog() {
  const input = document.getElementById("conceptSearchInput");
  if (!input) return;
  renderConceptSearch(input.value || "");
}

function selectWorkflowConceptFromButton(event) {
  const button = event.target.closest("[data-concept-id]");
  const concept = workflowConceptItems().find((candidate) => candidate.id === button?.dataset.conceptId);
  if (!concept) return;
  const query = document.getElementById("conceptSearchInput")?.value || concept.title || concept.id;
  const score = conceptSearchResults(query).find((result) => result.concept.id === concept.id)?.score || 0;
  renderConceptDocument(concept, score);
}

function renderConceptSearch(query) {
  const resultsEl = document.getElementById("conceptSearchResults");
  const documentEl = document.getElementById("conceptDocument");
  if (!resultsEl || !documentEl) return;
  const results = conceptSearchResults(query);
  if (!String(query || "").trim()) {
    resultsEl.innerHTML = "";
    documentEl.innerHTML = `<p class="concept-placeholder">Describe an analysis goal to find the closest workflow concept.</p>`;
    return;
  }
  if (!results.length) {
    resultsEl.innerHTML = `<p class="validation-bad">No workflow concept matched this description.</p>`;
    documentEl.innerHTML = `<p class="concept-placeholder">Try mentioning your inputs, desired outputs, or analysis type.</p>`;
    return;
  }
  resultsEl.innerHTML = results.slice(0, 4).map((result, index) => `
    <button class="concept-result ${index === 0 ? "active" : ""}" type="button" data-action="select-workflow-concept" data-concept-id="${escapeHtml(result.concept.id)}">
      <span>
        <strong>${escapeHtml(result.concept.title)}</strong>
        <small>${escapeHtml(result.concept.path || "")}</small>
      </span>
      <span class="concept-score">${Math.round(result.score * 100)}%</span>
    </button>
  `).join("");
  renderConceptDocument(results[0].concept, results[0].score);
}

function conceptSearchResults(query) {
  const concepts = workflowConceptItems();
  const queryVector = textVector(query);
  const queryMagnitude = vectorMagnitude(queryVector);
  if (!queryMagnitude) return [];
  const documentFrequencies = conceptDocumentFrequencies(concepts);
  return concepts
    .map((concept) => {
      const text = [
        concept.title,
        concept.id,
        concept.url,
        concept.path,
        concept.content
      ].join(" ");
      const docVector = textVector(text);
      const score = cosineSimilarity(queryVector, docVector, queryMagnitude, vectorMagnitude(docVector), documentFrequencies, concepts.length);
      return { concept, score };
    })
    .filter((result) => result.score > 0)
    .sort((a, b) => b.score - a.score || a.concept.title.localeCompare(b.concept.title));
}

function conceptDocumentFrequencies(concepts) {
  const frequencies = new Map();
  concepts.forEach((concept) => {
    const tokens = new Set(tokenizeConceptText([concept.title, concept.id, concept.content].join(" ")));
    tokens.forEach((token) => frequencies.set(token, (frequencies.get(token) || 0) + 1));
  });
  return frequencies;
}

function textVector(text) {
  return tokenizeConceptText(text).reduce((vector, token) => {
    vector.set(token, (vector.get(token) || 0) + 1);
    return vector;
  }, new Map());
}

function tokenizeConceptText(text) {
  const stopWords = new Set([
    "the", "and", "for", "with", "that", "this", "from", "into", "your", "you", "are", "can", "use", "uses",
    "workflow", "anvio", "anvi", "want", "need", "faire", "pour", "avec", "dans", "des", "les", "une", "mon", "mes",
    "and", "or", "to", "of", "in", "a", "is", "it", "on", "by", "as", "be", "when", "what", "how", "have", "has", "had"
  ]);
  return String(text || "")
    .toLowerCase()
    .normalize("NFKD")
    .replace(/[\u0300-\u036f]/g, "")
    .split(/[^a-z0-9]+/)
    .filter((token) => (token.length > 2 || token === "db" || token === "fa") && !stopWords.has(token))
    .flatMap((token) => expandConceptToken(stemConceptToken(token)));
}

function stemConceptToken(token) {
  return token
    .replace(/(ing|ers|ies|ied|ed|es|s)$/i, "")
    .replace(/genomic$/, "genom")
    .replace(/metagenomic$/, "metagenom");
}

function expandConceptToken(token) {
  const synonyms = {
    db: ["database"],
    database: ["db"],
    databases: ["db", "database"],
    contig: ["contigs"],
    contigs: ["contig"],
    fasta: ["fasta", "fna", "fa"],
    fna: ["fasta", "fa"],
    fastq: ["fastq", "reads", "read"],
    read: ["read", "reads", "fastq"],
    reads: ["read", "reads", "fastq"],
    bam: ["bam", "mapping", "alignment"],
    sam: ["sam", "mapping", "alignment"],
    profile: ["profile", "coverage"],
    variability: ["variation", "variant", "snv"],
    variation: ["variability", "variant", "snv"],
    variants: ["variability", "variant", "snv"],
    snv: ["variant", "variability"],
    assemblage: ["assembly", "assemble", "contig"],
    assembler: ["assemble", "assembly", "contig"],
    lecture: ["read", "reads", "fastq"],
    lectures: ["read", "reads", "fastq"],
    mapping: ["map", "mapping", "bowtie"],
    mapper: ["map", "mapping", "bowtie"],
    couverture: ["coverage", "profile"],
    abondance: ["abundance", "profile"],
    taxonomie: ["taxonomy"],
    fonction: ["function", "annotation"],
    annotation: ["annotation", "function"],
    genome: ["genome", "genomes"],
    genomes: ["genome", "genomes"],
    comparaison: ["compare", "comparative"],
    comparer: ["compare", "comparative"],
    pangenome: ["pangenome", "pangenomics", "pan"],
    arbre: ["tree", "phylogeny"],
    phylogenie: ["tree", "phylogeny", "phylogenomics"],
    phylogenomie: ["tree", "phylogeny", "phylogenomics"],
    telecharger: ["download", "sra"],
    telechargement: ["download", "sra"],
    trna: ["trna", "trnaseq"],
    arn: ["rna", "trna"],
    ecologie: ["ecology", "ecophylo"],
    distribution: ["distribution", "ecophylo"],
    hmm: ["hmm", "hmms"]
  };
  return [token, ...(synonyms[token] || [])].map(stemConceptToken);
}

function vectorMagnitude(vector) {
  return Math.sqrt([...vector.values()].reduce((sum, value) => sum + value * value, 0));
}

function cosineSimilarity(queryVector, docVector, queryMagnitude, docMagnitude, documentFrequencies, documentCount) {
  if (!queryMagnitude || !docMagnitude) return 0;
  let dot = 0;
  queryVector.forEach((queryValue, token) => {
    const docValue = docVector.get(token) || 0;
    if (!docValue) return;
    const idf = Math.log((1 + documentCount) / (1 + (documentFrequencies.get(token) || 0))) + 1;
    dot += queryValue * docValue * idf * idf;
  });
  return Math.min(1, dot / (queryMagnitude * docMagnitude));
}

function renderConceptDocument(concept, score = 0) {
  const documentEl = document.getElementById("conceptDocument");
  if (!documentEl) return;
  documentEl.innerHTML = `
    <div class="concept-document-meta">
      <span class="chip workflow">${escapeHtml(concept.id)}</span>
      <span class="chip relation">match ${Math.round(score * 100)}%</span>
      ${concept.url ? `<a href="${escapeHtml(concept.url)}" target="_blank" rel="noopener">Open anvi'o documentation</a>` : ""}
    </div>
    ${markdownToHtml(concept.content || "")}
  `;
}

function markdownToHtml(markdown) {
  const lines = String(markdown || "").split(/\r?\n/);
  const html = [];
  let paragraph = [];
  let listOpen = false;
  let codeOpen = false;
  const flushParagraph = () => {
    if (!paragraph.length) return;
    html.push(`<p>${inlineMarkdownToHtml(paragraph.join(" "))}</p>`);
    paragraph = [];
  };
  const closeList = () => {
    if (!listOpen) return;
    html.push("</ul>");
    listOpen = false;
  };
  lines.forEach((line) => {
    if (line.startsWith("```")) {
      flushParagraph();
      closeList();
      if (codeOpen) {
        html.push("</code></pre>");
        codeOpen = false;
      } else {
        html.push("<pre><code>");
        codeOpen = true;
      }
      return;
    }
    if (codeOpen) {
      html.push(escapeHtml(line));
      return;
    }
    if (!line.trim()) {
      flushParagraph();
      closeList();
      return;
    }
    const heading = line.match(/^(#{1,4})\s+(.+)$/);
    if (heading) {
      flushParagraph();
      closeList();
      const level = Math.min(heading[1].length + 1, 5);
      html.push(`<h${level}>${inlineMarkdownToHtml(heading[2])}</h${level}>`);
      return;
    }
    const listItem = line.match(/^-\s+(.+)$/);
    if (listItem) {
      flushParagraph();
      if (!listOpen) {
        html.push("<ul>");
        listOpen = true;
      }
      html.push(`<li>${inlineMarkdownToHtml(listItem[1])}</li>`);
      return;
    }
    paragraph.push(line.trim());
  });
  flushParagraph();
  closeList();
  if (codeOpen) html.push("</code></pre>");
  return html.join("\n");
}

function inlineMarkdownToHtml(text) {
  return escapeHtml(text)
    .replace(/`([^`]+)`/g, "<code>$1</code>")
    .replace(/\[([^\]]+)\]\(([^)]+)\)/g, `<a href="$2" target="_blank" rel="noopener">$1</a>`);
}

function addDocumentedPipelineFromButton(event) {
  const button = event.target.closest("[data-pipeline-id]");
  const pipelineId = button?.dataset.pipelineId || "";
  const pipeline = findDocumentedPipeline(pipelineId);
  if (!pipeline) return;
  const point = visibleSpawnPoint();
  withHistory(() => {
    if (pipeline.steps) addPreloadedWorkflow(pipeline, point);
    else spawnItem(pipeline, point.x, point.y, true, true, true);
  });
  if (typeof dom.modal.close === "function") dom.modal.close();
}

function addDocumentedPipelineToNewTabFromButton(event) {
  const button = event.target.closest("[data-pipeline-id]");
  const pipelineId = button?.dataset.pipelineId || "";
  const pipeline = findDocumentedPipeline(pipelineId);
  if (!pipeline) return;
  const cleanName = compactText((pipeline.title || "Workflow").replace(/\s+workflow$/i, ""), 30) || "Workflow";
  const point = { x: WORLD_WIDTH / 2, y: WORLD_HEIGHT / 2 };
  withHistory(() => {
    const tab = createTab(cleanName);
    state.tabs.push(tab);
    state.activeTabId = tab.id;
    state.selectedNodeId = "";
    state.selectedEdgeId = "";
    state.connectStartId = "";
    if (pipeline.steps) addPreloadedWorkflow(pipeline, point);
    else spawnItem(pipeline, point.x, point.y, true, false, false);
  });
  if (typeof dom.modal.close === "function") dom.modal.close();
  requestAnimationFrame(centerCanvasViewport);
}

function addPreloadedWorkflow(workflow, point) {
  const tab = currentTab();
  if (!tab) return;
  const previousSelectedNodeId = state.selectedNodeId;
  const previousSelectedEdgeId = state.selectedEdgeId;
  state.selectedNodeId = "";
  state.selectedEdgeId = "";

  const steps = expandPreloadedWorkflowSteps(workflow);
  const nodesByStep = new Map();
  const columns = Math.max(1, Math.ceil(Math.sqrt(steps.length)));
  steps.forEach((step, index) => {
    const item = workflowStepToItem(step, workflow);
    const col = index % columns;
    const row = Math.floor(index / columns);
    const node = createNodeFromItem(item, point.x + col * 330, point.y + row * 190);
    node.workflowSource = workflow.source || "";
    node.workflowRule = step.rule || step.id;
    tab.nodes.push(node);
    nodesByStep.set(step.id, node);
  });

  const initialInputs = new Map();
  const produced = new Map();
  steps.forEach((step) => {
    const consumer = nodesByStep.get(step.id);
    (step.inputs || []).forEach((artifact) => {
      const producer = produced.get(artifact);
      if (producer && producer.id !== consumer.id) {
        ensureEdge(producer.id, consumer.id, artifact);
        return;
      }
      if (!initialInputs.has(artifact)) initialInputs.set(artifact, []);
      initialInputs.get(artifact).push(consumer);
    });
    (step.outputs || []).forEach((artifact) => {
      produced.set(artifact, consumer);
    });
  });

  [...initialInputs.entries()].forEach(([artifact, consumers], index) => {
    const dataNode = createNodeFromItem(findDataBlock(artifact), point.x - 360, point.y + index * 135);
    dataNode.value = artifact;
    tab.nodes.push(dataNode);
    consumers.forEach((consumer) => ensureEdge(dataNode.id, consumer.id, artifact));
  });

  state.selectedNodeId = previousSelectedNodeId;
  state.selectedEdgeId = previousSelectedEdgeId;
  applyOptimizedPlacement(point);
  requestAnimationFrame(centerCanvasViewport);
}

function expandPreloadedWorkflowSteps(workflow, seen = new Set()) {
  if (!workflow || seen.has(workflow.id)) return [];
  seen.add(workflow.id);
  const allWorkflows = Array.isArray(window.ANVIO_PRELOADED_WORKFLOWS) ? window.ANVIO_PRELOADED_WORKFLOWS : [];
  const included = (workflow.includes || [])
    .flatMap((id) => expandPreloadedWorkflowSteps(allWorkflows.find((candidate) => candidate.id === id), seen));
  const ownSteps = (workflow.steps || []).map((step) => {
    if (!step.ref) return step;
    const [workflowId, stepId] = step.ref.split(":");
    const sourceWorkflow = allWorkflows.find((candidate) => candidate.id === workflowId);
    const sourceStep = sourceWorkflow?.steps?.find((candidate) => candidate.id === stepId);
    return sourceStep ? { ...sourceStep, id: step.id || `${workflowId}-${stepId}` } : step;
  }).filter((step) => step.command);
  return [...included, ...ownSteps].map((step, index) => ({
    ...step,
    id: `${workflow.id}:${step.id || step.rule || index}`
  }));
}

function workflowStepToItem(step, workflow) {
  const libraryItem = currentVersion().programs.find((item) => item.command === step.command || item.id === step.command);
  const externalFromAnvio = !libraryItem || !String(step.command || "").startsWith("anvi-");
  const parameters = (step.params || []).map((name) => ({
    name,
    value: "",
    defaultValue: null,
    required: false,
    enabled: externalFromAnvio,
    source: "github-snakefile",
    documentation: `Parameter used by ${step.rule || step.id} in ${workflow.source || "an anvi'o Snakefile"}.`
  }));
  return {
    ...(libraryItem ? clone(libraryItem) : {}),
    id: `workflow-step:${workflow.id}:${step.id}`,
    kind: "program",
    title: step.rule || step.command,
    command: step.command,
    description: step.description || libraryItem?.description || `Rule from ${workflow.title}.`,
    requires: clone(step.inputs || libraryItem?.requires || []),
    provides: clone(step.outputs || libraryItem?.provides || []),
    parameters: libraryItem ? mergeWorkflowParameters(clone(libraryItem.parameters || []), parameters) : parameters,
    examples: [`${step.command} ${parameters.map((param) => param.name).join(" ")}`.trim()],
    docsUrl: libraryItem?.docsUrl || "https://github.com/merenlab/anvio/tree/master/anvio/workflows",
    documentation: `${step.description || ""}\n\nSource: ${workflow.source || "merenlab/anvio"}\nRule: ${step.rule || step.id}`.trim(),
    externalFromAnvio,
    externalSource: workflow.source || "anvi'o Snakefile"
  };
}

function mergeWorkflowParameters(baseParameters, workflowParameters) {
  const byName = new Map(baseParameters.map((param) => [param.name, param]));
  workflowParameters.forEach((param) => {
    if (!byName.has(param.name)) byName.set(param.name, param);
  });
  return [...byName.values()];
}

function handleModalBackdropClick(event) {
  if (event.target !== dom.modal && event.target.closest("#modalContent")) return;
  if (typeof dom.modal.close === "function") dom.modal.close();
}

function searchableItemText(item) {
  return [
    item.command,
    item.id,
    item.title,
    item.artifactType,
    item.kind === "data" ? artifactKeywordText(item) : commandKeywordText(item.command || item.id),
    ...(item.aliases || []),
    ...(item.requires || []),
    ...(item.provides || []),
    item.description
  ].join(" ").toLowerCase();
}

function searchResultTitle(item) {
  return item.command || item.id || item.title;
}

function scoreField(field, query, words, weights) {
  const text = String(field || "").toLowerCase();
  if (!text) return 0;
  const tokens = text.split(/[^a-z0-9]+/).filter(Boolean);
  let score = 0;

  if (text === query) score += weights.exact || 0;
  else if (text.startsWith(query)) score += weights.prefix || 0;
  else if (text.includes(query)) score += weights.phrase || 0;

  words.forEach((word) => {
    if (tokens.includes(word)) score += weights.wordExact || 0;
    else if (text.includes(word)) score += weights.wordPartial || 0;
  });

  return score;
}

function scoreSearchItem(item, query, words) {
  const searchable = searchableItemText(item);
  const itemKeywords = new Set(tokenizeConceptText(searchable));
  const coverage = words.length
    ? words.filter((word) => itemKeywords.has(word) || searchable.includes(word)).length / words.length
    : 0;
  if (words.length && coverage < Math.min(1, words.length <= 2 ? 1 : 0.72)) return 0;

  const command = item.command || "";
  const commandBase = command.replace(/^anvi-/, "");
  let score = 1 + coverage * 240;
  words.forEach((word) => {
    if (itemKeywords.has(word)) score += item.kind === "data" ? 42 : 28;
  });
  if (commandBase === query) score += 1400;
  else if (commandBase.startsWith(`${query}-`)) score += 900;
  score += scoreField(command, query, words, {
    exact: 3000,
    prefix: 1600,
    phrase: 900,
    wordExact: 220,
    wordPartial: 140
  });
  score += scoreField(item.id, query, words, {
    exact: 2400,
    prefix: 1300,
    phrase: 760,
    wordExact: 170,
    wordPartial: 110
  });
  score += scoreField(item.title, query, words, {
    exact: 1200,
    prefix: 680,
    phrase: 420,
    wordExact: 90,
    wordPartial: 55
  });
  score += scoreField([item.artifactType, ...(item.requires || []), ...(item.provides || [])].join(" "), query, words, {
    exact: 700,
    prefix: 420,
    phrase: 280,
    wordExact: 70,
    wordPartial: 42
  });
  score += scoreField(item.description, query, words, {
    exact: 180,
    prefix: 90,
    phrase: 60,
    wordExact: 18,
    wordPartial: 9
  });
  if (item.kind === "data") {
    const dataText = artifactKeywordText(item).toLowerCase();
    score += scoreField(dataText, query, words, {
      exact: 2600,
      prefix: 1450,
      phrase: 860,
      wordExact: 210,
      wordPartial: 95
    });
    if (/\b(db|database)\b/.test(query) && /-db$/.test(item.id || "")) score += 360;
    if (/\b(file|path)\b/.test(query) && /-files?$/.test(item.id || "")) score += 220;
  }

  if (item.kind === "program" && command) score += 35;
  return score;
}

function estimatePillRows(items, width) {
  const values = items.length ? items : ["none"];
  const available = Math.max(72, width - 18);
  let rows = 1;
  let rowWidth = 0;

  values.forEach((item) => {
    const pillWidth = clamp(String(item || "").length * 6.4 + 18, 42, available);
    const nextWidth = rowWidth ? rowWidth + 4 + pillWidth : pillWidth;
    if (nextWidth > available && rowWidth) {
      rows += 1;
      rowWidth = pillWidth;
    } else {
      rowWidth = nextWidth;
    }
  });

  return rows;
}

function requiredNodeHeight(node, width = node.w) {
  const header = 35;
  const footer = 31;
  const bodyPadding = 16;
  const rowHeight = 22;

  if (node.kind === "data") {
    const labels = [
      node.artifactType,
      node.value ? compactText(node.value, 42) : ""
    ].filter(Boolean);
    return Math.ceil(header + footer + bodyPadding + estimatePillRows(labels, width) * rowHeight);
  }

  const inputRows = estimatePillRows(node.requires || [], width);
  const outputRows = estimatePillRows(node.provides || [], width);
  const modifiedRows = modifiedParameterLabels(node).length
    ? estimatePillRows(modifiedParameterLabels(node), width)
    : 0;
  const groupLabelHeight = 13;
  const groupGap = 7;
  const groupInternalGap = 4;
  const groupCount = 2 + (modifiedRows ? 1 : 0);
  const bodyHeight = bodyPadding
    + groupLabelHeight * groupCount
    + groupInternalGap * groupCount
    + groupGap * Math.max(0, groupCount - 1)
    + (inputRows + outputRows + modifiedRows) * rowHeight;

  return Math.ceil(header + footer + bodyHeight);
}

function normalizeNodeSize(node) {
  const minimumHeight = requiredNodeHeight(node);
  node.h = Math.max(node.h || 0, minimumHeight, node.kind === "data" ? DATA_NODE_DEFAULT_HEIGHT : NODE_DEFAULT_HEIGHT);
  node.y = clamp(node.y, 0, WORLD_HEIGHT - node.h);
}

function normalizeActiveNodeSizes() {
  activeNodes().forEach(normalizeNodeSize);
}

function renderNodes() {
  normalizeActiveNodeSizes();
  const nodes = activeNodes();
  dom.nodeLayer.innerHTML = nodes.map((node) => `
    <div class="block-node ${node.kind} ${node.externalFromAnvio ? "external-tool" : ""} ${node.id === state.selectedNodeId ? "selected" : ""}"
      data-node-id="${node.id}"
      style="left:${node.x}px; top:${node.y}px; width:${node.w}px; height:${node.h}px">
      <div class="block-head">
        <span class="block-kind"></span>
        <span class="block-title" title="${escapeHtml(node.title)}">${escapeHtml(nodeDisplayTitle(node))}</span>
      </div>
      <div class="block-body">${renderNodeArtifacts(node)}</div>
      <div class="block-foot">
        <span>${escapeHtml(node.kind === "data" ? node.title : node.kind)}</span>
      </div>
      <span class="port port-top" data-port-side="top" title="Connect from top"></span>
      <span class="port port-right" data-port-side="right" title="Connect from right"></span>
      <span class="port port-bottom" data-port-side="bottom" title="Connect from bottom"></span>
      <span class="port port-left" data-port-side="left" title="Connect from left"></span>
      <div class="resize-handle" data-resize-node="${node.id}" title="Resize"></div>
    </div>
  `).join("");
}

function isExternalAnvioItem(item) {
  return Boolean(item?.externalFromAnvio || (item?.command && !String(item.command).startsWith("anvi-") && item?.source === "github-snakefile"));
}

function nodeDisplayTitle(node) {
  if (node.kind === "data") return node.refId || node.title;
  return node.command || node.title;
}

function renderNodeArtifacts(node) {
  if (node.kind === "data") {
    const labels = [
      { type: "data", text: node.refId || node.title },
      node.artifactType ? { type: "data", text: node.artifactType, artifact: node.refId || node.title } : null,
      node.value ? { type: "value", text: compactText(node.value, 42) } : null
    ].filter(Boolean);
    return labels.length
      ? `<div class="artifact-strip">${labels.map((label) => renderArtifactPill(label.type, label.text, label.artifact)).join("")}</div>`
      : `<div class="artifact-empty">Data artifact</div>`;
  }

  const inputs = node.requires || [];
  const outputs = node.provides || [];
  const modifiedParameters = modifiedParameterLabels(node);
  return `
    <div class="artifact-group input-group">
      <span class="artifact-label">Inputs</span>
      <div class="artifact-strip">
        ${inputs.length ? inputs.map((artifact) => renderInputArtifactPill(node, artifact)).join("") : `<span class="artifact-empty">none</span>`}
      </div>
    </div>
    <div class="artifact-group output-group">
      <span class="artifact-label">Outputs</span>
      <div class="artifact-strip">
        ${outputs.length ? outputs.map((artifact) => renderArtifactPill("output", artifact)).join("") : `<span class="artifact-empty">none</span>`}
      </div>
    </div>
    ${modifiedParameters.length ? `
      <div class="artifact-group params-group">
        <span class="artifact-label">Modified params</span>
        <div class="artifact-strip">
          ${modifiedParameters.map((parameter) => renderArtifactPill("parameter", parameter)).join("")}
        </div>
      </div>
    ` : ""}
  `;
}

function renderInputArtifactPill(node, artifact) {
  const state = inputArtifactCompatibilityState(node, artifact);
  return renderArtifactPill("input", artifact, artifact, state.status, state.reason);
}

function renderArtifactPill(type, text, artifact = text, status = "", reason = "") {
  const isArtifact = ["input", "output", "data"].includes(type);
  if (isArtifact && artifact) {
    const statusClass = status ? ` ${status}` : "";
    const statusText = reason ? `${reason} ` : "";
    return `<button class="artifact-pill ${type}${statusClass}" type="button" data-artifact-label="${escapeHtml(artifact)}" data-artifact-role="${escapeHtml(type === "data" ? "both" : type)}" title="${escapeHtml(statusText)}Show related blocks for ${escapeHtml(artifact)}">${escapeHtml(text)}</button>`;
  }
  return `<span class="artifact-pill ${type}" title="${escapeHtml(text)}">${escapeHtml(text)}</span>`;
}

function handleNodeLayerClick(event) {
  const pill = event.target.closest("[data-artifact-label]");
  if (!pill) return;
  event.preventDefault();
  event.stopPropagation();
  const artifact = pill.dataset.artifactLabel;
  const role = pill.dataset.artifactRole || "both";
  if (!artifact) return;
  state.focusedArtifact = {
    artifact,
    role,
    nodeId: pill.closest("[data-node-id]")?.dataset.nodeId || ""
  };
  if (dom.searchInput) dom.searchInput.value = "";
  state.searchQuery = "";
  renderSearch();
}

function modifiedParameterLabels(node) {
  if (!node || node.kind === "data") return [];
  const inputLabels = Object.values(node.inputValues || {})
    .filter((input) => String(input.value || "").trim())
    .map((input) => `${input.flag || `--${input.artifact}`}=${compactText(input.value, 36)}`);

  const parameterLabels = (node.parameters || [])
    .filter((param) => parameterIsModified(param))
    .map((param) => {
      const value = String(param.value ?? "").trim();
      return value ? `${param.name}=${compactText(value, 36)}` : param.name;
    });

  return [...inputLabels, ...parameterLabels];
}

function parameterIsModified(param) {
  if (!param) return false;
  const value = String(param.value ?? "").trim();
  const defaultValue = param.defaultValue == null ? "" : String(param.defaultValue).trim();
  if (param.source === "user") return true;
  if (param.enabled && !param.required) return true;
  if (value && value !== defaultValue) return true;
  return false;
}


function renderEdges() {
  const tab = currentTab();
  if (!tab) return;
  const pieces = [`
    <defs>
      <marker id="arrowhead" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto" markerUnits="strokeWidth">
        <path d="M0,0 L10,4 L0,8 z"></path>
      </marker>
    </defs>
  `];

  tab.edges.forEach((edge) => {
    const points = routeEdge(edge);
    if (points.length < 2) return;
    const d = pointsToPath(points);
    const isSelected = edge.id === state.selectedEdgeId;
    const isHovered = edge.id === state.hoveredEdgeId;
    const className = `edge-path ${isSelected ? "selected" : ""} ${isHovered ? "hovered" : ""}`.trim();
    pieces.push(`<g class="edge-group" data-edge-id="${edge.id}">`);
    pieces.push(`<path class="${className}" data-edge-id="${edge.id}" d="${d}"></path>`);
    pieces.push(`<path class="edge-hit" data-edge-id="${edge.id}" d="${d}"></path>`);

    if (isSelected) {
      const middle = points[Math.floor(points.length / 2)];
      const x = middle.x;
      const y = middle.y - 38;
      pieces.push(`<g class="edge-toolbar" data-edge-id="${edge.id}">
        <rect class="edge-toolbar-bg" x="${x - 48}" y="${y - 16}" width="96" height="32" rx="7"></rect>
        <text class="edge-toolbar-btn" data-edge-action="add-point" x="${x - 28}" y="${y + 1}" text-anchor="middle" dominant-baseline="middle">+</text>
        <text class="edge-toolbar-btn" data-edge-action="reset-edge" x="${x}" y="${y + 1}" text-anchor="middle" dominant-baseline="middle">R</text>
        <text class="edge-toolbar-btn danger" data-edge-action="delete-edge" x="${x + 28}" y="${y + 1}" text-anchor="middle" dominant-baseline="middle">x</text>
      </g>`);
      points.slice(1, -1).forEach((point, index) => {
        pieces.push(`<circle class="edge-handle" data-edge-id="${edge.id}" data-point-index="${index}" cx="${point.x}" cy="${point.y}" r="5"></circle>`);
      });
    }
    pieces.push("</g>");
  });

  dom.edgeLayer.innerHTML = pieces.join("");
}

function renderToolbarState() {
  document.querySelectorAll('[data-action="undo"]').forEach((button) => {
    button.disabled = state.history.past.length === 0;
  });
  document.querySelectorAll('[data-action="redo"]').forEach((button) => {
    button.disabled = state.history.future.length === 0;
  });
  document.querySelectorAll("[data-theme-switch]").forEach((switchInput) => {
    switchInput.checked = Boolean(state.settings.nightMode);
  });
  document.querySelectorAll('[data-action="toggle-library-panel"]').forEach((button) => {
    const collapsed = Boolean(state.settings.libraryCollapsed);
    button.setAttribute("aria-expanded", String(!collapsed));
    button.title = collapsed ? "Expand search panel" : "Collapse search panel";
    button.setAttribute("aria-label", button.title);
  });
  document.querySelectorAll('[data-action="toggle-inspector-panel"]').forEach((button) => {
    const collapsed = Boolean(state.settings.inspectorCollapsed);
    button.setAttribute("aria-expanded", String(!collapsed));
    button.title = collapsed ? "Expand inspector panel" : "Collapse inspector panel";
    button.setAttribute("aria-label", button.title);
  });
  if (dom.zoomLabel) dom.zoomLabel.textContent = `${Math.round(canvasZoom() * 100)}%`;
}

function renderInspector() {
  const node = activeNodes().find((candidate) => candidate.id === state.selectedNodeId);
  const edge = activeEdges().find((candidate) => candidate.id === state.selectedEdgeId);
  if (node) {
    renderNodeInspector(node);
    return;
  }
  if (edge) {
    renderEdgeInspector(edge);
    return;
  }
  dom.inspector.innerHTML = "";
  if (!state.settings.inspectorCollapsed) {
    state.settings.inspectorCollapsed = true;
    applyPanelState();
    localStorage.setItem("anvio-builder-settings", JSON.stringify(state.settings));
    renderToolbarState();
  }
}

function renderNodeInspector(node) {
  enforceParameterCompatibility(node);
  const libraryItem = findLibraryItem(node.kind, node.refId);
  const title = libraryItem?.title || node.title;
  const description = libraryItem?.description || node.description || "No description available.";
  const documentation = libraryItem?.documentation || node.documentation || description;
  const requiredParameters = (node.parameters || []).filter((param) => param.required);
  const optionalParameters = (node.parameters || []).filter((param) => !param.required);
  const compatibility = parameterCompatibilityState(node);
  ensureInputValues(node);

  dom.inspector.innerHTML = `
    <div class="inspector-content">
      <div class="inspector-section readonly-summary">
        <h2>${escapeHtml(title)}</h2>
        <p class="doc-description">${escapeHtml(description)}</p>
        <div class="chips">
          <span class="chip ${node.kind}">${escapeHtml(node.kind)}</span>
          ${node.externalFromAnvio ? `<span class="chip external">external from anvi'o</span>` : ""}
          ${libraryItem?.artifactType ? `<span class="chip data">${escapeHtml(libraryItem.artifactType)}</span>` : ""}
          <span class="chip">${escapeHtml(currentVersion().label)}</span>
        </div>
        ${node.externalFromAnvio ? `<p class="external-note">External command from ${escapeHtml(node.externalSource || "a workflow")}. Parameters are editable and custom flags can be added below.</p>` : ""}
      </div>

      ${node.kind === "data" ? `
        <div class="inspector-section">
          <h3>Data path</h3>
          <label class="field-stack">
            <span class="small-label">Path for ${escapeHtml(node.refId || node.title || "dataset")}</span>
            <input data-editable data-node-field="value" value="${escapeHtml(node.value || "")}" placeholder="${escapeHtml(artifactPathPlaceholder(node.refId || node.title, node.artifactType))}">
            <small class="param-doc">Use a file path, Galaxy dataset path, or workflow variable that resolves to this data block.</small>
          </label>
        </div>
      ` : renderRequiredInputSection(node, requiredParameters, compatibility)}

      ${node.kind !== "data" ? `
        <div class="inspector-section advanced-section">
          <button class="advanced-toggle" type="button" data-action="toggle-advanced">
            <span>${node.externalFromAnvio ? "External and custom parameters" : "Advanced parameters"}</span>
            <span class="toggle-arrow">▶</span>
          </button>
          <div class="advanced-content" ${node.externalFromAnvio ? "" : "hidden"}>
            ${optionalParameters.length ? `
              <div class="parameters-list">
                ${optionalParameters.map((param) => {
                  const globalIndex = node.parameters.indexOf(param);
                  const lock = compatibility.locks.get(globalIndex);
                  return `
                    <div class="param-row ${lock ? "locked" : ""}">
                      <input data-param-enabled="${globalIndex}" type="checkbox" ${param.enabled ? "checked" : ""} ${lock ? "disabled" : ""} aria-label="Enable ${escapeHtml(param.name)}">
                      <label>
                        <input class="param-name-input" data-editable data-param-name="${globalIndex}" value="${escapeHtml(param.name)}" placeholder="--parameter" ${lock ? "disabled" : ""}>
                        <input data-editable data-param-value="${globalIndex}" value="${escapeHtml(param.value || "")}" placeholder="${param.defaultValue ? `default: ${escapeHtml(param.defaultValue)}` : "value"}" ${lock ? "disabled" : ""}>
                        <small>${escapeHtml(cleanDocText(param.documentation) || param.source || "Custom parameter")}</small>
                        ${lock ? `<small class="param-lock">${escapeHtml(lock.reason)}</small>` : ""}
                      </label>
                    </div>
                  `;
                }).join("")}
              </div>
            ` : `<p>No additional parameters are documented for this command.</p>`}
            <div class="inspector-actions">
              <button class="secondary-button" type="button" data-action="add-param">+ Add parameter</button>
              <button class="secondary-button" type="button" data-action="add-variable">+ Add variable</button>
            </div>
            ${renderVariableEditor(Object.entries(node.variables || {}))}
          </div>
        </div>
      ` : ""}

      <div class="inspector-section">
        <h3>Complete documentation</h3>
        <pre class="full-documentation">${escapeHtml(documentation)}</pre>
      </div>

      <div class="inspector-section">
        <h3>Official documentation</h3>
        <button class="doc-button" type="button" data-action="open-official-doc" data-doc-url="${escapeHtml(libraryItem?.docsUrl || node.docsUrl || currentVersion().docsUrl)}">Show official anvi'o documentation</button>
      </div>

      <div class="inspector-actions">
        <button class="danger-icon-button" type="button" data-action="delete-node" title="Delete block" aria-label="Delete block">
          <svg viewBox="0 0 24 24" aria-hidden="true"><path d="M3 6h18"></path><path d="M8 6V4h8v2"></path><path d="M19 6l-1 15H6L5 6"></path><path d="M10 11v6"></path><path d="M14 11v6"></path></svg>
        </button>
      </div>
    </div>
  `;
}

function renderRequiredInputSection(node, parameters, compatibility = parameterCompatibilityState(node)) {
  const inputs = Object.values(ensureInputValues(node));
  return `
    <div class="inspector-section">
      <h3>Documented inputs and required parameters</h3>
      ${inputs.length ? `
        <div class="parameters-list">
          ${inputs.map((input) => `
            <label class="field-stack param-input-row">
              <span class="param-header">
                <code class="param-name">${escapeHtml(input.flag)}</code>
                <span class="param-default">${escapeHtml(input.artifact)}</span>
              </span>
              <input data-editable data-input-value="${escapeHtml(input.artifact)}" value="${escapeHtml(input.value || "")}" placeholder="${escapeHtml(artifactPathPlaceholder(input.artifact))}">
              <small class="param-doc">Input path or dataset variable for ${escapeHtml(input.artifact)}.</small>
            </label>
          `).join("")}
        </div>
      ` : ""}
      ${parameters.length ? `
        <div class="parameters-list">
          ${parameters.map((param) => {
            const index = (node.parameters || []).indexOf(param);
            const defaultValue = param.defaultValue ?? "";
            const value = param.value || defaultValue || "";
            const lock = compatibility.locks.get(index);
            return `
              <label class="field-stack param-input-row ${lock ? "locked" : ""}">
                <span class="param-header">
                  <code class="param-name">${escapeHtml(param.name)}</code>
                  ${defaultValue !== "" ? `<span class="param-default">default: ${escapeHtml(defaultValue)}</span>` : `<span class="param-default">required</span>`}
                </span>
                <input data-editable data-required-param-value="${index}" value="${escapeHtml(value)}" placeholder="Enter value for ${escapeHtml(param.name)}" ${lock ? "disabled" : ""}>
                ${cleanDocText(param.documentation) || param.source ? `<small class="param-doc">${escapeHtml(cleanDocText(param.documentation) || param.source)}</small>` : ""}
                ${lock ? `<small class="param-lock">${escapeHtml(lock.reason)}</small>` : ""}
              </label>
            `;
          }).join("")}
        </div>
      ` : ""}
      ${!inputs.length && !parameters.length ? `<p>No command inputs are listed in the local runtime database for this item.</p>` : ""}
      <div class="command-preview"><code id="commandPreview">${escapeHtml(buildCommand(node))}</code></div>
    </div>
  `;
}

function renderParameterEditor(parameters) {
  return `
    <h3>Parameters</h3>
    <div class="inspector-section">
      ${parameters.length ? parameters.map((param, index) => `
        <div class="param-row">
          <input data-param-enabled="${index}" type="checkbox" ${param.enabled ? "checked" : ""} aria-label="Enable ${escapeHtml(param.name)}">
          <label>
            <input class="param-name-input" data-editable data-param-name="${index}" value="${escapeHtml(param.name)}" placeholder="--parameter">
            <input data-editable data-param-value="${index}" value="${escapeHtml(param.value || "")}" placeholder="${param.defaultValue ? `default: ${escapeHtml(param.defaultValue)}` : "value"}">
            <small>${escapeHtml(cleanDocText(param.documentation) || param.source || "Custom parameter")}</small>
          </label>
        </div>
      `).join("") : `<p>No parameters were explicitly visible in the official help page for this command. Add custom parameters as needed.</p>`}
    </div>
  `;
}

function renderVariableEditor(variables) {
  return `
    <h3>Variables</h3>
    <div class="inspector-section">
      ${variables.length ? variables.map(([name, value], index) => `
        <div class="var-row">
          <button class="danger-button" type="button" data-variable-remove="${escapeHtml(name)}" aria-label="Remove variable">-</button>
          <label>
            <input data-editable data-variable-name="${index}" value="${escapeHtml(name)}" placeholder="name">
            <input data-editable data-variable-value="${escapeHtml(name)}" value="${escapeHtml(value)}" placeholder="value">
          </label>
        </div>
      `).join("") : `<p>Variables are exported with the workspace and can be used to remember Galaxy-friendly dataset names, paths, or sample ids.</p>`}
    </div>
  `;
}

function renderEdgeInspector(edge) {
  const nodes = activeNodes();
  const fromNode = nodeById(edge.from);
  const toNode = nodeById(edge.to);
  dom.inspector.innerHTML = `
    <h2>Arrow</h2>
    <p>${escapeHtml(fromNode?.title || "Unknown")} -> ${escapeHtml(toNode?.title || "Unknown")}</p>
    <div class="inspector-section">
      <label class="field-stack">
        <span class="small-label">From</span>
        <select data-edge-field="from">
          ${nodes.map((node) => `<option value="${node.id}" ${node.id === edge.from ? "selected" : ""}>${escapeHtml(node.title)}</option>`).join("")}
        </select>
      </label>
      <label class="field-stack">
        <span class="small-label">To</span>
        <select data-edge-field="to">
          ${nodes.map((node) => `<option value="${node.id}" ${node.id === edge.to ? "selected" : ""}>${escapeHtml(node.title)}</option>`).join("")}
        </select>
      </label>
      <label class="field-stack">
        <span class="small-label">Label</span>
        <input data-editable data-edge-field="label" value="${escapeHtml(edge.label || "")}" placeholder="optional artifact or note">
      </label>
    </div>
    <p>Select the arrow to show its canvas toolbar. Use + to add an inflexion point, R to reset the path, x to delete the arrow, and drag red handles to reshape the path. Shift-click a handle to remove that point.</p>
    <div class="inspector-actions">
      <button class="secondary-button" type="button" data-action="reset-edge">Reset path</button>
      <button class="danger-button" type="button" data-action="delete-edge">Delete arrow</button>
    </div>
  `;
}

function handleInspectorInput(event) {
  const node = activeNodes().find((candidate) => candidate.id === state.selectedNodeId);
  const edge = activeEdges().find((candidate) => candidate.id === state.selectedEdgeId);
  const target = event.target;

  if (node && target.dataset.nodeField) {
    node[target.dataset.nodeField] = target.value;
    if (target.dataset.nodeField === "command") updateCommandPreview(node);
    if (["title", "description", "command", "value"].includes(target.dataset.nodeField)) {
      normalizeNodeSize(node);
      renderNodes();
      renderEdges();
    }
  }
  if (node && target.dataset.nodeNumber) {
    node[target.dataset.nodeNumber] = Number(target.value);
    renderNodes();
    renderEdges();
  }
  if (node && target.dataset.requiredParamValue) {
    const param = node.parameters[Number(target.dataset.requiredParamValue)];
    if (param) {
      param.value = target.value;
      param.enabled = true;
      enforceParameterCompatibility(node, param.name);
      updateCommandPreview(node);
      normalizeNodeSize(node);
      renderNodes();
      renderEdges();
    }
  }
  if (node && target.dataset.inputValue) {
    ensureInputValues(node);
    const input = node.inputValues[target.dataset.inputValue];
    if (input) {
      input.value = target.value;
      enforceParameterCompatibility(node);
      updateCommandPreview(node);
      normalizeNodeSize(node);
      renderNodes();
      renderEdges();
    }
  }
  if (node && target.dataset.paramName) {
    const param = node.parameters[Number(target.dataset.paramName)];
    if (param) {
      param.name = target.value;
      enforceParameterCompatibility(node, param.name);
      updateCommandPreview(node);
      normalizeNodeSize(node);
      renderNodes();
      renderEdges();
    }
  }
  if (node && target.dataset.paramValue) {
    const param = node.parameters[Number(target.dataset.paramValue)];
    if (param) {
      param.value = target.value;
      enforceParameterCompatibility(node, param.name);
      updateCommandPreview(node);
      normalizeNodeSize(node);
      renderNodes();
      renderEdges();
    }
  }
  if (node && target.dataset.variableValue) {
    node.variables[target.dataset.variableValue] = target.value;
  }
  if (node && target.dataset.variableName) {
    const entries = Object.entries(node.variables);
    const index = Number(target.dataset.variableName);
    const oldName = entries[index]?.[0];
    if (oldName && target.value && target.value !== oldName) {
      const value = node.variables[oldName];
      delete node.variables[oldName];
      node.variables[target.value] = value;
    }
  }
  if (edge && target.dataset.edgeField === "label") {
    edge.label = target.value;
  }
}

function handleInspectorChange(event) {
  const node = activeNodes().find((candidate) => candidate.id === state.selectedNodeId);
  const edge = activeEdges().find((candidate) => candidate.id === state.selectedEdgeId);
  const target = event.target;

  if (node && target.dataset.paramEnabled) {
    const param = node.parameters[Number(target.dataset.paramEnabled)];
    if (param) {
      pushHistory();
      param.enabled = target.checked;
      enforceParameterCompatibility(node, param.name);
      updateCommandPreview(node);
      normalizeNodeSize(node);
      render();
    }
  }
  if (edge && target.dataset.edgeField && target.dataset.edgeField !== "label") {
    pushHistory();
    edge[target.dataset.edgeField] = target.value;
    edge.points = [];
    delete edge.layoutRoute;
    render();
  }
  if (target.dataset.variableRemove) {
    withHistory(() => {
      delete node.variables[target.dataset.variableRemove];
    });
  }
}

function updateCommandPreview(node) {
  const preview = document.getElementById("commandPreview");
  if (preview) preview.textContent = buildCommand(node);
}

function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

function canvasZoom() {
  const zoom = Number(state.settings.zoom);
  return Number.isFinite(zoom) ? clamp(zoom, MIN_ZOOM, MAX_ZOOM) : 1;
}

function applyCanvasZoom() {
  const zoom = canvasZoom();
  state.settings.zoom = zoom;
  dom.canvasWorld.style.width = `${Math.round(WORLD_WIDTH * zoom)}px`;
  dom.canvasWorld.style.height = `${Math.round(WORLD_HEIGHT * zoom)}px`;
  dom.edgeLayer.style.transform = `scale(${zoom})`;
  dom.nodeLayer.style.transform = `scale(${zoom})`;
  if (dom.zoomLabel) dom.zoomLabel.textContent = `${Math.round(zoom * 100)}%`;
}

function viewportWorldCenter() {
  const zoom = canvasZoom();
  return {
    x: (dom.canvasFrame.scrollLeft + dom.canvasFrame.clientWidth / 2) / zoom,
    y: (dom.canvasFrame.scrollTop + dom.canvasFrame.clientHeight / 2) / zoom
  };
}

function setCanvasZoom(value, center = viewportWorldCenter()) {
  const nextZoom = clamp(Number(value.toFixed(2)), MIN_ZOOM, MAX_ZOOM);
  if (nextZoom === canvasZoom()) return;
  state.settings.zoom = nextZoom;
  applySettings();
  dom.canvasFrame.scrollLeft = Math.max(0, Math.round(center.x * nextZoom - dom.canvasFrame.clientWidth / 2));
  dom.canvasFrame.scrollTop = Math.max(0, Math.round(center.y * nextZoom - dom.canvasFrame.clientHeight / 2));
  renderToolbarState();
}

function zoomCanvas(delta, center = viewportWorldCenter()) {
  setCanvasZoom(canvasZoom() + delta, center);
}

function centerCanvasViewport() {
  const zoom = canvasZoom();
  const frame = dom.canvasFrame;
  frame.scrollLeft = Math.max(0, Math.round((WORLD_WIDTH * zoom - frame.clientWidth) / 2));
  frame.scrollTop = Math.max(0, Math.round((WORLD_HEIGHT * zoom - frame.clientHeight) / 2));
}

function visibleSpawnPoint() {
  const zoom = canvasZoom();
  const x = clamp(
    (dom.canvasFrame.scrollLeft + dom.canvasFrame.clientWidth / 2) / zoom - NODE_DEFAULT_WIDTH / 2,
    20,
    WORLD_WIDTH - 360
  );
  const y = clamp(
    (dom.canvasFrame.scrollTop + dom.canvasFrame.clientHeight / 2) / zoom - NODE_DEFAULT_HEIGHT / 2,
    20,
    WORLD_HEIGHT - 220
  );
  return { x, y };
}

function laneY(originY, index, total, spacing = 172) {
  return Math.round(originY + (index - (Math.max(1, total) - 1) / 2) * spacing);
}

function rectsOverlap(a, b, margin = 26) {
  return !(
    a.x + a.w + margin < b.x ||
    b.x + b.w + margin < a.x ||
    a.y + a.h + margin < b.y ||
    b.y + b.h + margin < a.y
  );
}

function findOpenNodePosition(x, y, w = NODE_DEFAULT_WIDTH, h = NODE_DEFAULT_HEIGHT, ignoreIds = new Set()) {
  const base = {
    x: clamp(Math.round(x), 20, WORLD_WIDTH - w - 20),
    y: clamp(Math.round(y), 20, WORLD_HEIGHT - h - 20)
  };
  const existing = activeNodes().filter((node) => !ignoreIds.has(node.id));
  const offsets = [{ x: 0, y: 0 }];
  const stepX = 292;
  const stepY = 168;

  for (let ring = 1; ring <= 6; ring += 1) {
    for (let dx = -ring; dx <= ring; dx += 1) {
      for (let dy = -ring; dy <= ring; dy += 1) {
        if (Math.abs(dx) !== ring && Math.abs(dy) !== ring) continue;
        offsets.push({ x: dx * stepX, y: dy * stepY });
      }
    }
  }

  const candidates = offsets
    .map((offset) => ({
      x: clamp(base.x + offset.x, 20, WORLD_WIDTH - w - 20),
      y: clamp(base.y + offset.y, 20, WORLD_HEIGHT - h - 20),
      distance: Math.hypot(offset.x, offset.y)
    }))
    .sort((a, b) => a.distance - b.distance || a.y - b.y || a.x - b.x);

  return candidates.find((candidate) => {
    const rect = { x: candidate.x, y: candidate.y, w, h };
    return !existing.some((node) => rectsOverlap(rect, node));
  }) || base;
}

function normalizedDataKey(value) {
  return String(value || "")
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, "-")
    .replace(/^-+|-+$/g, "");
}

function dataNodeMatchesItem(node, item) {
  if (node.kind !== "data") return false;
  if (node.refId && item.id && node.refId === item.id) return true;
  if (node.artifactType && item.artifactType && node.artifactType === item.artifactType) {
    const nodeKeys = new Set([
      normalizedDataKey(node.refId),
      normalizedDataKey(node.title)
    ]);
    return [
      normalizedDataKey(item.id),
      normalizedDataKey(item.title)
    ].some((key) => key && nodeKeys.has(key));
  }
  return false;
}

function connectDirectFunctionBlocks(producer, consumer) {
  if (!producer || !consumer || producer.id === consumer.id) return;
  if (producer.kind === "data" || consumer.kind === "data") return;

  (consumer.requires || []).forEach((requiredArtifact) => {
    if ((producer.provides || []).some((providedArtifact) => artifactCompatible(requiredArtifact, providedArtifact))) {
      ensureEdge(producer.id, consumer.id, requiredArtifact);
    }
  });
}

function connectExistingFunctionBlocks(node) {
  if (!node || node.kind === "data") return;
  activeNodes()
    .filter((candidate) => candidate.id !== node.id && candidate.kind !== "data")
    .forEach((candidate) => {
      connectDirectFunctionBlocks(candidate, node);
      connectDirectFunctionBlocks(node, candidate);
    });
}

function connectDataBlockToExistingConsumers(dataNode) {
  if (!dataNode || dataNode.kind !== "data" || !dataNode.refId) return;
  activeNodes()
    .filter((candidate) => candidate.id !== dataNode.id && candidate.kind !== "data")
    .filter((candidate) => (candidate.requires || []).some((artifact) => artifactCompatible(artifact, dataNode.refId)))
    .forEach((consumer) => ensureEdge(dataNode.id, consumer.id, dataNode.refId));
}

function inferredEdgeLabel(fromNode, toNode) {
  if (!fromNode || !toNode || fromNode.id === toNode.id) return "";
  if (fromNode.kind === "data" && toNode.kind !== "data") {
    return (toNode.requires || []).find((artifact) => artifactCompatible(artifact, fromNode.refId)) || "";
  }
  if (fromNode.kind !== "data" && toNode.kind === "data") {
    return (fromNode.provides || []).find((artifact) => artifactCompatible(toNode.refId, artifact)) || "";
  }
  if (fromNode.kind !== "data" && toNode.kind !== "data") {
    return (toNode.requires || []).find((required) => (fromNode.provides || []).some((provided) => artifactCompatible(required, provided))) || "";
  }
  return "";
}

function ensureInferredEdge(fromNode, toNode) {
  const label = inferredEdgeLabel(fromNode, toNode);
  if (!label) return null;
  return ensureEdge(fromNode.id, toNode.id, label);
}

function ensureManualEdge(fromNode, toNode, options = {}) {
  if (!fromNode || !toNode || fromNode.id === toNode.id) return null;
  const edge = ensureInferredEdge(fromNode, toNode) || ensureEdge(fromNode.id, toNode.id, "");
  if (edge && options.fromSide) edge.fromSide = options.fromSide;
  if (edge && options.toSide) edge.toSide = options.toSide;
  return edge;
}

function dataNodeHasProducer(node) {
  return activeEdges().some((edge) => edge.to === node.id && nodeById(edge.from)?.kind !== "data");
}

function producerForDataNode(node) {
  const edge = activeEdges().find((candidate) => candidate.to === node.id && nodeById(candidate.from)?.kind !== "data");
  return edge ? nodeById(edge.from) : null;
}

function reusableDataNode(item, x, y) {
  const matches = activeNodes().filter((node) => dataNodeMatchesItem(node, item));
  if (!matches.length) return null;
  return matches
    .map((node) => ({
      node,
      produced: dataNodeHasProducer(node) ? 0 : 1,
      distance: Math.hypot(node.x - x, node.y - y)
    }))
    .sort((a, b) => a.produced - b.produced || a.distance - b.distance)[0].node;
}

function nodeProvidesArtifact(node, artifactId) {
  if (!node) return false;
  if (node.kind === "data") return artifactCompatible(artifactId, node.refId);
  return (node.provides || []).some((providedArtifact) => artifactCompatible(artifactId, providedArtifact));
}

function nodeRequiresArtifact(node, artifactId) {
  return node?.kind !== "data" && (node.requires || []).some((requiredArtifact) => artifactCompatible(requiredArtifact, artifactId));
}

function existingProducerForArtifact(artifactId, consumerId = "") {
  return activeNodes()
    .filter((node) => node.id !== consumerId && node.kind !== "data" && nodeProvidesArtifact(node, artifactId))
    .sort((a, b) => a.x - b.x || a.y - b.y)[0] || null;
}

function connectRequiredArtifact(node, artifactId, index, total) {
  const existingProducer = existingProducerForArtifact(artifactId, node.id);
  if (existingProducer) {
    ensureEdge(existingProducer.id, node.id, artifactId);
    return;
  }

  const y = laneY(node.y, index, total);
  const dataItem = findDataBlock(artifactId);
  if (!dataItem) return;

  const existingDataNode = reusableDataNode(dataItem, node.x - 330, y);
  if (existingDataNode) {
    const producerNode = producerForDataNode(existingDataNode);
    ensureEdge((producerNode || existingDataNode).id, node.id, artifactId);
  }
}

function connectProvidedArtifact(node, artifactId) {
  activeNodes()
    .filter((candidate) => candidate.id !== node.id && nodeRequiresArtifact(candidate, artifactId))
    .forEach((consumer) => ensureEdge(node.id, consumer.id, inferredEdgeLabel(node, consumer) || artifactId));

  const dataItem = findDataBlock(artifactId);
  if (!dataItem) return;
  const existingDataNode = reusableDataNode(dataItem, node.x + 330, node.y);
  if (existingDataNode && !dataNodeHasProducer(existingDataNode)) {
    ensureEdge(node.id, existingDataNode.id, artifactId);
  }
}

function spawnItem(item, x, y, withIO = true, selectCreated = true, autoLayout = true) {
  if (!item) return null;

  const tab = currentTab();
  if (!tab) {
    showModal("Error", "No active canvas. Please create a new canvas first.");
    return null;
  }

  if (tab.versionId && item.versionId && tab.versionId !== item.versionId) {
    showModal("Version Mismatch", `
      <p>This block is from <strong>${item.versionLabel || item.versionId}</strong>, but the current canvas is locked to <strong>${getVersionLabel(tab.versionId)}</strong>.</p>
      <p>Cannot add blocks from different versions to the same canvas.</p>
    `);
    return null;
  }

  const node = createNodeFromItem(item, x, y);
  tab.nodes.push(node);

  if (withIO && item.kind !== "data") {
    const inputIds = item.requires || [];
    const outputIds = item.provides || [];
    inputIds.forEach((artifactId, index) => {
      connectRequiredArtifact(node, artifactId, index, inputIds.length);
    });
    outputIds.forEach((artifactId) => {
      connectProvidedArtifact(node, artifactId);
    });
    connectExistingFunctionBlocks(node);
  } else if (withIO && item.kind === "data") {
    connectDataBlockToExistingConsumers(node);
  }
  if (selectCreated) {
    state.selectedNodeId = node.id;
    state.selectedEdgeId = "";
    state.focusedArtifact = null;
    if (dom.searchInput) dom.searchInput.value = "";
    state.searchQuery = "";
    revealInspectorPanel();
  }
  if (autoLayout && activeNodes().length > 1) {
    applyOptimizedPlacement(activeGraphCenter({ x, y }));
  }
  return node;
}

function createNodeFromItem(item, x, y) {
  const baseHeight = item.kind === "data" ? DATA_NODE_DEFAULT_HEIGHT : NODE_DEFAULT_HEIGHT;
  const h = Math.max(baseHeight, requiredNodeHeight({
    kind: item.kind,
    artifactType: item.artifactType || "",
    value: "",
    requires: item.requires || [],
    provides: item.provides || [],
    w: NODE_DEFAULT_WIDTH
  }, NODE_DEFAULT_WIDTH));
  const position = findOpenNodePosition(x, y, NODE_DEFAULT_WIDTH, h);
  const node = {
    id: uid("node"),
    refId: item.id,
    kind: item.kind,
    title: item.title || item.id,
    command: item.command || "",
    artifactType: item.artifactType || "",
    description: item.description || "",
    docsUrl: item.docsUrl || officialDocsUrlForItem(item.kind === "data" ? "artifact" : item.kind, item.command || item.id, currentVersion().upstreamVersionId || currentVersion().id || "main"),
    requires: clone(item.requires || []),
    provides: clone(item.provides || []),
    parameters: clone(item.parameters || []).map((param) => ({
      ...param,
      enabled: Boolean(param.enabled || param.required || item.externalFromAnvio),
      value: param.value || (param.required && param.defaultValue != null ? String(param.defaultValue) : "")
    })),
    examples: clone(item.examples || []),
    documentation: item.documentation || "",
    externalFromAnvio: Boolean(item.externalFromAnvio),
    externalSource: item.externalSource || "",
    inputValues: {},
    variables: {},
    value: "",
    x: position.x,
    y: position.y,
    w: NODE_DEFAULT_WIDTH,
    h
  };
  ensureInputValues(node);
  normalizeNodeSize(node);
  return node;
}

function ensureEdge(from, to, label = "") {
  const existing = activeEdges().find((edge) => edge.from === from && edge.to === to && edge.label === label);
  if (existing) return existing;
  if (label) {
    const unlabeled = activeEdges().find((edge) => edge.from === from && edge.to === to && !edge.label);
    if (unlabeled) {
      unlabeled.label = label;
      return unlabeled;
    }
  }
  const edge = {
    id: uid("edge"),
    from,
    to,
    label,
    points: []
  };
  currentTab().edges.push(edge);
  return edge;
}

function nodeById(id) {
  return activeNodes().find((node) => node.id === id);
}

function handleNodePointerDown(event) {
  if (event.target.closest("[data-artifact-label]")) return;
  const resize = event.target.closest("[data-resize-node]");
  const port = event.target.closest(".port");
  const nodeElement = event.target.closest("[data-node-id]");
  if (!nodeElement) return;
  const node = nodeById(nodeElement.dataset.nodeId);
  if (!node) return;

  if (state.mode === "connect") {
    event.preventDefault();
    handleConnectNode(node);
    return;
  }

  // Handle arrow creation from port drag
  if (port) {
    event.preventDefault();
    event.stopPropagation();
    pushHistory();
    nodeElement.setPointerCapture(event.pointerId);

    const startPoint = canvasPoint(event);
    let previewLine = null;
    let targetNode = null;
    let targetPort = null;

    const onMove = (moveEvent) => {
      const currentPoint = canvasPoint(moveEvent);

      // Create or update preview line
      if (!previewLine) {
        previewLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        previewLine.setAttribute('x1', startPoint.x);
        previewLine.setAttribute('y1', startPoint.y);
        previewLine.setAttribute('stroke', '#145f68');
        previewLine.setAttribute('stroke-width', '2');
        previewLine.setAttribute('stroke-dasharray', '5,5');
        previewLine.setAttribute('class', 'connection-preview');
        previewLine.setAttribute('pointer-events', 'none');
        dom.edgeLayer.appendChild(previewLine);
      }
      previewLine.setAttribute('x2', currentPoint.x);
      previewLine.setAttribute('y2', currentPoint.y);

      // Find closest port within reasonable distance
      const allPorts = Array.from(document.querySelectorAll('.port'));
      let closest = { port: null, node: null, distance: 100 };

      allPorts.forEach(p => {
        const blockNode = p.closest('[data-node-id]');
        if (!blockNode || blockNode.dataset.nodeId === node.id) return;

        const rect = p.getBoundingClientRect();
        const portX = rect.left + rect.width / 2;
        const portY = rect.top + rect.height / 2;
        const distance = Math.hypot(moveEvent.clientX - portX, moveEvent.clientY - portY);

        if (distance < closest.distance) {
          closest = { port: p, node: blockNode, distance };
        }
      });

      let closestBlock = { node: null, distance: 80 };
      document.querySelectorAll("[data-node-id]").forEach((blockNode) => {
        if (blockNode.dataset.nodeId === node.id) return;
        const rect = blockNode.getBoundingClientRect();
        const dx = Math.max(rect.left - moveEvent.clientX, 0, moveEvent.clientX - rect.right);
        const dy = Math.max(rect.top - moveEvent.clientY, 0, moveEvent.clientY - rect.bottom);
        const distance = Math.hypot(dx, dy);
        if (distance < closestBlock.distance) closestBlock = { node: blockNode, distance };
      });

      if (closest.port) {
        targetNode = nodeById(closest.node.dataset.nodeId);
        targetPort = closest.port;
        closest.port.style.background = '#145f68';
      } else if (closestBlock.node) {
        if (targetPort) targetPort.style.background = '';
        targetNode = nodeById(closestBlock.node.dataset.nodeId);
        targetPort = null;
      } else if (targetPort) {
        targetPort.style.background = '';
        targetNode = null;
        targetPort = null;
      } else {
        targetNode = null;
      }
    };

    const onUp = () => {
      nodeElement.releasePointerCapture(event.pointerId);
      nodeElement.removeEventListener("pointermove", onMove);
      nodeElement.removeEventListener("pointerup", onUp);

      // Clean up preview line
      previewLine?.remove();
      if (targetPort) targetPort.style.background = '';

      // Create arrow if valid target
      if (targetNode) {
        withHistory(() => {
          const edge = ensureManualEdge(node, targetNode, {
            fromSide: port.dataset.portSide || "",
            toSide: targetPort?.dataset.portSide || ""
          });
          state.selectedEdgeId = edge?.id || "";
          render();
        });
      }
      renderInspector();
    };

    nodeElement.addEventListener("pointermove", onMove);
    nodeElement.addEventListener("pointerup", onUp);
    return;
  }

  state.selectedNodeId = node.id;
  state.selectedEdgeId = "";
  state.focusedArtifact = null;
  if (dom.searchInput) dom.searchInput.value = "";
  state.searchQuery = "";
  revealInspectorPanel();
  renderSearch();
  renderInspector();
  renderEdges();

  const start = canvasPoint(event);
  const origin = { x: node.x, y: node.y, w: node.w, h: node.h };
  pushHistory();
  nodeElement.setPointerCapture(event.pointerId);

  const onMove = (moveEvent) => {
    const point = canvasPoint(moveEvent);
    if (resize) {
      node.w = Math.max(150, Math.min(520, origin.w + point.x - start.x));
      node.h = Math.max(requiredNodeHeight(node), Math.min(900, origin.h + point.y - start.y));
    } else {
      node.x = Math.max(0, Math.min(WORLD_WIDTH - node.w, origin.x + point.x - start.x));
      node.y = Math.max(0, Math.min(WORLD_HEIGHT - node.h, origin.y + point.y - start.y));
    }
    nodeElement.style.left = `${node.x}px`;
    nodeElement.style.top = `${node.y}px`;
    nodeElement.style.width = `${node.w}px`;
    nodeElement.style.height = `${node.h}px`;
    renderEdges();
  };
  const onUp = () => {
    nodeElement.releasePointerCapture(event.pointerId);
    nodeElement.removeEventListener("pointermove", onMove);
    nodeElement.removeEventListener("pointerup", onUp);
    renderInspector();
  };
  nodeElement.addEventListener("pointermove", onMove);
  nodeElement.addEventListener("pointerup", onUp);
}

function handleConnectNode(node) {
  if (!state.connectStartId) {
    state.connectStartId = node.id;
    state.selectedNodeId = node.id;
    if (dom.searchInput) dom.searchInput.value = "";
    state.searchQuery = "";
    revealInspectorPanel();
    renderSearch();
    renderInspector();
    renderNodes();
    return;
  }
  if (state.connectStartId === node.id) {
    state.connectStartId = "";
    return;
  }
  withHistory(() => {
    const startNode = nodeById(state.connectStartId);
    const edge = ensureManualEdge(startNode, node);
    state.selectedNodeId = "";
    state.selectedEdgeId = edge?.id || "";
    state.mode = "select";
    state.connectStartId = "";
  });
}


function handleEdgePointerDown(event) {
  const toolbarAction = event.target.closest("[data-edge-action]");
  const handle = event.target.closest(".edge-handle");
  const edgeId = event.target.closest("[data-edge-id]")?.dataset.edgeId;
  if (!edgeId) return;
  const edge = activeEdges().find((candidate) => candidate.id === edgeId);
  if (!edge) return;

  event.preventDefault();
  event.stopPropagation();
  state.selectedEdgeId = edge.id;
  state.selectedNodeId = "";

  if (toolbarAction) {
    const action = toolbarAction.dataset.edgeAction;

    if (action === "add-point") {
      withHistory(() => {
        addEdgePoint(edge, canvasPoint(event));
      });
      return;
    }

    if (action === "reset-edge") {
      withHistory(() => {
        edge.points = [];
        delete edge.layoutRoute;
      });
      return;
    }

    if (action === "delete-edge") {
      withHistory(() => {
        const tab = currentTab();
        tab.edges = tab.edges.filter((candidate) => candidate.id !== edgeId);
        state.selectedEdgeId = "";
      });
      return;
    }
  }

  if (handle) {
    const pointIndex = Number(handle.dataset.pointIndex);
    if (!Number.isFinite(pointIndex)) return;
    if (!edge.points.length) edge.points = routeEdge(edge).slice(1, -1);
    if (event.shiftKey) {
      withHistory(() => {
        edge.points.splice(pointIndex, 1);
      });
      return;
    }
    pushHistory();
    const pointerId = event.pointerId;

    const onMove = (moveEvent) => {
      if (moveEvent.pointerId !== pointerId) return;
      moveEvent.preventDefault();
      edge.points[pointIndex] = canvasPoint(moveEvent);
      renderEdges();
    };
    const onUp = (upEvent) => {
      if (upEvent.pointerId !== pointerId) return;
      document.removeEventListener("pointermove", onMove);
      document.removeEventListener("pointerup", onUp);
      document.removeEventListener("pointercancel", onUp);
      renderInspector();
    };
    document.addEventListener("pointermove", onMove);
    document.addEventListener("pointerup", onUp);
    document.addEventListener("pointercancel", onUp);
    return;
  }

  if (event.altKey) {
    withHistory(() => {
      addEdgePoint(edge, canvasPoint(event));
    });
    return;
  }

  render();
}

function handleCanvasWorldPointerDown(event) {
  if (event.button !== 0) return;
  const target = event.target;
  if (!(target instanceof Element)) return;
  if (target.closest(".intent-builder")) return;
  if (target.closest("[data-node-id], [data-edge-id], [data-edge-action], .edge-handle")) return;
  if (!state.selectedNodeId && !state.selectedEdgeId) return;
  state.selectedNodeId = "";
  state.selectedEdgeId = "";
  render();
}

function handleOutsideSelectionPointerDown(event) {
  if (event.button !== 0 || !state.selectedEdgeId) return;
  const target = event.target;
  if (!(target instanceof Element)) return;
  if (dom.canvasWorld.contains(target) || dom.inspectorPanel.contains(target) || dom.modal.contains(target)) return;
  state.selectedEdgeId = "";
  render();
}

function handleCanvasWheel(event) {
  if (!event.ctrlKey && !event.metaKey) return;
  event.preventDefault();
  const zoom = canvasZoom();
  const rect = dom.canvasWorld.getBoundingClientRect();
  const center = {
    x: (event.clientX - rect.left) / zoom,
    y: (event.clientY - rect.top) / zoom
  };
  zoomCanvas(event.deltaY < 0 ? ZOOM_IN_STEP : -ZOOM_OUT_STEP, center);
}

function handleCanvasPanPointerDown(event) {
  if (event.button !== 0) return;
  const target = event.target;
  if (!(target instanceof Element)) return;
  if (target.closest(".intent-builder")) return;
  if (target.closest("[data-node-id], [data-edge-id], [data-edge-action], .edge-handle")) return;

  event.preventDefault();
  const pointerId = event.pointerId;
  const start = {
    x: event.clientX,
    y: event.clientY,
    scrollLeft: dom.canvasFrame.scrollLeft,
    scrollTop: dom.canvasFrame.scrollTop
  };
  dom.canvasFrame.classList.add("panning");

  const onMove = (moveEvent) => {
    if (moveEvent.pointerId !== pointerId) return;
    moveEvent.preventDefault();
    dom.canvasFrame.scrollLeft = start.scrollLeft - (moveEvent.clientX - start.x);
    dom.canvasFrame.scrollTop = start.scrollTop - (moveEvent.clientY - start.y);
  };
  const onUp = (upEvent) => {
    if (upEvent.pointerId !== pointerId) return;
    dom.canvasFrame.classList.remove("panning");
    document.removeEventListener("pointermove", onMove);
    document.removeEventListener("pointerup", onUp);
    document.removeEventListener("pointercancel", onUp);
  };
  document.addEventListener("pointermove", onMove);
  document.addEventListener("pointerup", onUp);
  document.addEventListener("pointercancel", onUp);
}

function addEdgePoint(edge, point) {
  const route = routeEdge(edge);
  const interior = edge.points.length ? [...edge.points] : route.slice(1, -1);
  const insertAt = nearestSegmentIndex(route, point);
  interior.splice(Math.max(0, insertAt - 1), 0, point);
  edge.points = interior;
  delete edge.layoutRoute;
}

function nearestSegmentIndex(points, point) {
  let bestIndex = 1;
  let bestDistance = Infinity;
  for (let i = 0; i < points.length - 1; i += 1) {
    const distance = distanceToSegment(point, points[i], points[i + 1]);
    if (distance < bestDistance) {
      bestDistance = distance;
      bestIndex = i + 1;
    }
  }
  return bestIndex;
}

function distanceToSegment(point, a, b) {
  const dx = b.x - a.x;
  const dy = b.y - a.y;
  if (dx === 0 && dy === 0) return Math.hypot(point.x - a.x, point.y - a.y);
  const t = Math.max(0, Math.min(1, ((point.x - a.x) * dx + (point.y - a.y) * dy) / (dx * dx + dy * dy)));
  return Math.hypot(point.x - (a.x + t * dx), point.y - (a.y + t * dy));
}

function canvasPoint(event) {
  const rect = dom.canvasWorld.getBoundingClientRect();
  const zoom = canvasZoom();
  return {
    x: Math.round((event.clientX - rect.left) / zoom),
    y: Math.round((event.clientY - rect.top) / zoom)
  };
}

function routeEdge(edge) {
  const from = nodeById(edge.from);
  const to = nodeById(edge.to);
  if (!from || !to) return [];
  const source = edge.fromSide ? nodeSidePortInfo(from, edge.fromSide) : nodePortInfoToward(from, nodeCenter(to));
  const target = edge.toSide ? nodeSidePortInfo(to, edge.toSide) : nodePortInfoToward(to, nodeCenter(from));
  const defaultStart = source.point;
  const defaultEnd = target.point;
  if (edge.points?.length) {
    const start = edge.fromSide ? defaultStart : nodeBorderPointToward(from, edge.points[0] || nodeCenter(to));
    const end = edge.toSide ? defaultEnd : nodeBorderPointToward(to, edge.points[edge.points.length - 1] || nodeCenter(from));
    return simplifyRoutePoints([start, ...edge.points, end]);
  }
  return dynamicOrthogonalRoute(from, to, edge) || perpendicularFallbackRoute(source, target, from.id, to.id);
}

function attachLayoutRouteToNodes(route, from, to) {
  const interior = route.slice(1, -1).map((point) => ({ x: point.x, y: point.y }));
  const sourceTarget = interior[0] || route[1] || nodeCenter(to);
  const targetSource = interior[interior.length - 1] || route[route.length - 2] || nodeCenter(from);
  const start = nodePortTowardPoint(from, sourceTarget);
  const end = nodePortTowardPoint(to, targetSource);
  const attached = [start];
  interior.forEach((point) => pfgoAppendOrthogonalPoint(attached, point));
  pfgoAppendOrthogonalPoint(attached, end);
  return simplifyRoutePoints(attached);
}

function nodeCenter(node) {
  return {
    x: node.x + node.w / 2,
    y: node.y + node.h / 2
  };
}

function nodePortTowardPoint(node, point) {
  return nodeBorderPointToward(node, point);
}

function nodeBorderPointToward(node, point) {
  return nodePortInfoToward(node, point).point;
}

function nodeSideCenter(node, side) {
  const center = nodeCenter(node);
  if (side === "left") return { x: Math.round(node.x), y: Math.round(center.y) };
  if (side === "right") return { x: Math.round(node.x + node.w), y: Math.round(center.y) };
  if (side === "top") return { x: Math.round(center.x), y: Math.round(node.y) };
  return { x: Math.round(center.x), y: Math.round(node.y + node.h) };
}

function nodeSideNormal(side) {
  if (side === "left") return { x: -1, y: 0 };
  if (side === "right") return { x: 1, y: 0 };
  if (side === "top") return { x: 0, y: -1 };
  return { x: 0, y: 1 };
}

function nodeSidePortInfo(node, side, penalty = 0) {
  return {
    side,
    point: nodeSideCenter(node, side),
    normal: nodeSideNormal(side),
    penalty
  };
}

function nodePortInfoToward(node, point) {
  const center = nodeCenter(node);
  const dx = point.x - center.x;
  const dy = point.y - center.y;
  if (dx === 0 && dy === 0) {
    return nodeSidePortInfo(node, "right");
  }
  const side = Math.abs(dx) >= Math.abs(dy)
    ? (dx >= 0 ? "right" : "left")
    : (dy >= 0 ? "bottom" : "top");
  return nodeSidePortInfo(node, side);
}

function dynamicOrthogonalRoute(from, to, edge = {}) {
  const targetCenter = nodeCenter(to);
  const source = edge.fromSide ? nodeSidePortInfo(from, edge.fromSide) : nodePortInfoToward(from, targetCenter);
  const targetPorts = edge.toSide
    ? [nodeSidePortInfo(to, edge.toSide)]
    : nodeDestinationPorts(to, source.point);
  const obstacles = activeNodes().filter((node) => node.id !== from.id && node.id !== to.id);
  let best = null;

  targetPorts.forEach((target) => {
    const sourceStub = offsetPoint(source.point, source.normal, 36);
    const targetStub = offsetPoint(target.point, target.normal, 36);
    candidateOrthogonalRoutes(sourceStub, targetStub, obstacles).forEach((middleRoute) => {
      const route = [source.point, sourceStub, ...middleRoute.slice(1, -1), targetStub, target.point];
      const simplified = simplifyRoutePoints(route);
      const score = routeLength(simplified) + routeTurnCount(simplified) * 42 + target.penalty;
      if (!best || score < best.score) best = { route: simplified, score };
    });
  });

  return best?.route || null;
}

function perpendicularFallbackRoute(source, target, fromId, toId) {
  const sourceStub = offsetPoint(source.point, source.normal, 36);
  const targetStub = offsetPoint(target.point, target.normal, 36);
  const candidates = candidateOrthogonalRoutes(
    sourceStub,
    targetStub,
    activeNodes().filter((node) => node.id !== fromId && node.id !== toId)
  );
  const middle = candidates[0] || [sourceStub, { x: targetStub.x, y: sourceStub.y }, targetStub];
  return simplifyRoutePoints([source.point, sourceStub, ...middle.slice(1, -1), targetStub, target.point]);
}

function nodeDestinationPorts(node, fromPoint) {
  const preferred = nodePortInfoToward(node, fromPoint);
  const ports = ["left", "right", "top", "bottom"].map((side) => nodeSidePortInfo(node, side));
  const seen = new Set();
  return ports
    .map((port) => {
      const point = nodeSideCenter(node, port.side);
      return {
        ...port,
        point,
        penalty: (port.penalty || 0) + Math.hypot(point.x - preferred.point.x, point.y - preferred.point.y) * 0.8
      };
    })
    .filter((port) => {
      const key = routePointKey(port.point);
      if (seen.has(key)) return false;
      seen.add(key);
      return true;
    })
    .sort((a, b) => a.penalty - b.penalty || a.side.localeCompare(b.side));
}

function offsetPoint(point, normal, distance) {
  return {
    x: Math.round(point.x + normal.x * distance),
    y: Math.round(point.y + normal.y * distance)
  };
}

function candidateOrthogonalRoutes(start, end, obstacles) {
  const rawRoutes = [];
  const add = (points) => {
    const route = simplifyRoutePoints(points);
    if (route.length < 2) return;
    if (routeEverySegmentAllowed(route, obstacles)) rawRoutes.push(route);
  };

  add([start, { x: end.x, y: start.y }, end]);
  add([start, { x: start.x, y: end.y }, end]);

  const laneXs = routingLaneValues(start.x, end.x, obstacles, "x");
  const laneYs = routingLaneValues(start.y, end.y, obstacles, "y");
  laneXs.forEach((x) => {
    add([start, { x, y: start.y }, { x, y: end.y }, end]);
  });
  laneYs.forEach((y) => {
    add([start, { x: start.x, y }, { x: end.x, y }, end]);
  });

  laneXs.slice(0, 8).forEach((x) => {
    laneYs.slice(0, 8).forEach((y) => {
      add([start, { x, y: start.y }, { x, y }, { x: end.x, y }, end]);
      add([start, { x: start.x, y }, { x, y }, { x, y: end.y }, end]);
    });
  });

  return rawRoutes.sort((a, b) => routeLength(a) - routeLength(b) || routeTurnCount(a) - routeTurnCount(b)).slice(0, 24);
}

function routingLaneValues(startValue, endValue, obstacles, axis) {
  const limit = axis === "x" ? WORLD_WIDTH : WORLD_HEIGHT;
  const pad = 34;
  const values = new Set([
    Math.round((startValue + endValue) / 2),
    Math.round(clamp(startValue - 80, 20, limit - 20)),
    Math.round(clamp(startValue + 80, 20, limit - 20)),
    Math.round(clamp(endValue - 80, 20, limit - 20)),
    Math.round(clamp(endValue + 80, 20, limit - 20))
  ]);
  obstacles.forEach((node) => {
    const min = axis === "x" ? node.x : node.y;
    const max = axis === "x" ? node.x + node.w : node.y + node.h;
    [min - pad, max + pad, (min + max) / 2].forEach((value) => {
      values.add(Math.round(clamp(value, 20, limit - 20)));
    });
  });
  return [...values]
    .sort((a, b) => Math.abs(a - (startValue + endValue) / 2) - Math.abs(b - (startValue + endValue) / 2) || a - b)
    .slice(0, 12);
}

function routeEverySegmentAllowed(route, obstacles) {
  for (let index = 0; index < route.length - 1; index += 1) {
    if (!routeSegmentAllowed(route[index], route[index + 1], obstacles)) return false;
  }
  return true;
}

function routeBetweenPorts(start, end, obstacles) {
  const xs = routingAxisValues(start.x, end.x, obstacles, "x");
  const ys = routingAxisValues(start.y, end.y, obstacles, "y");
  const startKey = routePointKey(start);
  const endKey = routePointKey(end);
  const queue = [start];
  const parents = new Map([[startKey, null]]);
  const distance = new Map([[startKey, 0]]);

  while (queue.length) {
    const point = queue.shift();
    const key = routePointKey(point);
    if (key === endKey) break;
    routingNeighbors(point, xs, ys).forEach((next) => {
      const nextKey = routePointKey(next);
      if (parents.has(nextKey)) return;
      if (!routeSegmentAllowed(point, next, obstacles)) return;
      parents.set(nextKey, point);
      distance.set(nextKey, (distance.get(key) || 0) + Math.abs(next.x - point.x) + Math.abs(next.y - point.y));
      queue.push(next);
    });
    queue.sort((a, b) => {
      const da = (distance.get(routePointKey(a)) || 0) + Math.abs(a.x - end.x) + Math.abs(a.y - end.y);
      const db = (distance.get(routePointKey(b)) || 0) + Math.abs(b.x - end.x) + Math.abs(b.y - end.y);
      return da - db;
    });
  }

  if (!parents.has(endKey)) return null;
  const route = [];
  let cursor = end;
  while (cursor) {
    route.unshift(cursor);
    cursor = parents.get(routePointKey(cursor));
  }
  return route;
}

function routingAxisValues(startValue, endValue, obstacles, axis) {
  const values = new Set([
    Math.round(startValue),
    Math.round(endValue),
    axis === "x" ? 20 : 20,
    axis === "x" ? WORLD_WIDTH - 20 : WORLD_HEIGHT - 20,
    Math.round((startValue + endValue) / 2)
  ]);
  const pad = 28;
  obstacles.forEach((node) => {
    if (axis === "x") {
      [node.x - pad, node.x, node.x + node.w, node.x + node.w + pad, node.x + node.w / 2].forEach((value) => values.add(Math.round(clamp(value, 20, WORLD_WIDTH - 20))));
    } else {
      [node.y - pad, node.y, node.y + node.h, node.y + node.h + pad, node.y + node.h / 2].forEach((value) => values.add(Math.round(clamp(value, 20, WORLD_HEIGHT - 20))));
    }
  });
  [startValue - 72, startValue + 72, endValue - 72, endValue + 72].forEach((value) => {
    values.add(Math.round(clamp(value, 20, axis === "x" ? WORLD_WIDTH - 20 : WORLD_HEIGHT - 20)));
  });
  return [...values].sort((a, b) => a - b);
}

function routingNeighbors(point, xs, ys) {
  const xi = xs.indexOf(Math.round(point.x));
  const yi = ys.indexOf(Math.round(point.y));
  const neighbors = [];
  if (xi > 0) neighbors.push({ x: xs[xi - 1], y: point.y });
  if (xi >= 0 && xi < xs.length - 1) neighbors.push({ x: xs[xi + 1], y: point.y });
  if (yi > 0) neighbors.push({ x: point.x, y: ys[yi - 1] });
  if (yi >= 0 && yi < ys.length - 1) neighbors.push({ x: point.x, y: ys[yi + 1] });
  return neighbors;
}

function routeSegmentAllowed(a, b, obstacles) {
  if (a.x !== b.x && a.y !== b.y) return false;
  return !obstacles.some((node) => segmentIntersectsNode(a, b, node));
}

function routePointKey(point) {
  return `${Math.round(point.x)},${Math.round(point.y)}`;
}

function routeLength(points) {
  let total = 0;
  for (let index = 1; index < points.length; index += 1) {
    total += Math.abs(points[index].x - points[index - 1].x) + Math.abs(points[index].y - points[index - 1].y);
  }
  return total;
}

function routeTurnCount(points) {
  let turns = 0;
  for (let index = 2; index < points.length; index += 1) {
    const a = { x: Math.sign(points[index - 1].x - points[index - 2].x), y: Math.sign(points[index - 1].y - points[index - 2].y) };
    const b = { x: Math.sign(points[index].x - points[index - 1].x), y: Math.sign(points[index].y - points[index - 1].y) };
    if (a.x !== b.x || a.y !== b.y) turns += 1;
  }
  return turns;
}

function orthogonalRoute(start, end, fromId, toId) {
  const midX = Math.round((start.x + end.x) / 2);
  const direct = [
    start,
    { x: midX, y: start.y },
    { x: midX, y: end.y },
    end
  ];
  if (!routeIntersectsBlocks(direct, fromId, toId)) return direct;

  const obstacles = activeNodes().filter((node) => node.id !== fromId && node.id !== toId);
  const topLimit = Math.max(24, Math.min(start.y, end.y) - 120);
  const bottomLimit = Math.min(WORLD_HEIGHT - 24, Math.max(start.y, end.y) + 120);
  const topBlocked = obstacles.some((node) => topLimit > node.y - 24 && topLimit < node.y + node.h + 24);
  const detourY = topBlocked ? bottomLimit : topLimit;
  const startX = Math.min(WORLD_WIDTH - 20, start.x + 44);
  const endX = Math.max(20, end.x - 44);
  return [
    start,
    { x: startX, y: start.y },
    { x: startX, y: detourY },
    { x: endX, y: detourY },
    { x: endX, y: end.y },
    end
  ];
}

function routeIntersectsBlocks(points, fromId, toId) {
  const obstacles = activeNodes().filter((node) => node.id !== fromId && node.id !== toId);
  for (let i = 0; i < points.length - 1; i += 1) {
    const a = points[i];
    const b = points[i + 1];
    if (obstacles.some((node) => segmentIntersectsNode(a, b, node))) return true;
  }
  return false;
}

function segmentIntersectsNode(a, b, node) {
  const pad = 16;
  const left = node.x - pad;
  const right = node.x + node.w + pad;
  const top = node.y - pad;
  const bottom = node.y + node.h + pad;
  if (a.x === b.x) {
    const y1 = Math.min(a.y, b.y);
    const y2 = Math.max(a.y, b.y);
    return a.x >= left && a.x <= right && y2 >= top && y1 <= bottom;
  }
  if (a.y === b.y) {
    const x1 = Math.min(a.x, b.x);
    const x2 = Math.max(a.x, b.x);
    return a.y >= top && a.y <= bottom && x2 >= left && x1 <= right;
  }
  return false;
}

function pointsToPath(points) {
  return points.map((point, index) => `${index === 0 ? "M" : "L"} ${Math.round(point.x)} ${Math.round(point.y)}`).join(" ");
}

function buildCommand(node) {
  if (node.kind === "data") return `${node.refId || node.title}${node.value ? `=${node.value}` : ""}`;
  const parts = [node.command || node.title];
  selectedCommandInputs(node).forEach((input) => {
    parts.push(input.flag || `--${input.artifact}`);
    parts.push(shellQuote(input.resolvedValue));
  });
  (node.parameters || []).forEach((param) => {
    if (!param.enabled) return;
    parts.push(param.name);
    const value = param.value || (param.required && param.defaultValue != null ? String(param.defaultValue) : "");
    if (value) parts.push(shellQuote(value));
  });
  return parts.join(" ");
}

function shellQuote(value) {
  const text = String(value);
  if (!text) return "";
  if (/^[A-Za-z0-9_./:=+-]+$/.test(text)) return text;
  return `'${text.replaceAll("'", "'\"'\"'")}'`;
}

function deleteSelectedNode() {
  const nodeId = state.selectedNodeId;
  if (!nodeId) return;
  withHistory(() => {
    const tab = currentTab();
    tab.nodes = tab.nodes.filter((node) => node.id !== nodeId);
    state.selectedNodeId = "";
  });
}

function deleteSelectedEdge() {
  const edgeId = state.selectedEdgeId;
  if (!edgeId) return;
  withHistory(() => {
    const tab = currentTab();
    tab.edges = tab.edges.filter((edge) => edge.id !== edgeId);
    state.selectedEdgeId = "";
  });
}

function resetSelectedEdge() {
  const edge = activeEdges().find((candidate) => candidate.id === state.selectedEdgeId);
  if (!edge) return;
  withHistory(() => {
    edge.points = [];
    delete edge.layoutRoute;
  });
}

function addCustomParam() {
  const node = activeNodes().find((candidate) => candidate.id === state.selectedNodeId);
  if (!node || node.kind === "data") return;
  withHistory(() => {
    node.parameters.push({
      name: "--custom-parameter",
      value: "",
      defaultValue: null,
      required: false,
      enabled: true,
      source: "user",
      documentation: "User-added parameter."
    });
  });
}

function addVariable() {
  const node = activeNodes().find((candidate) => candidate.id === state.selectedNodeId);
  if (!node) return;
  withHistory(() => {
    let index = Object.keys(node.variables).length + 1;
    while (`variable_${index}` in node.variables) index += 1;
    node.variables[`variable_${index}`] = "";
  });
}

function validateWorkflow() {
  return validateWorkflowReport().flatMap((entry) => entry.issues.map((issue) => `${entry.title}: ${issue}`));
}

function validateWorkflowReport() {
  const report = activeNodes().map((node) => ({
    node,
    title: node.title || node.refId || "Untitled block",
    kind: node.kind || "program",
    issues: []
  }));
  const byId = new Map(report.map((entry) => [entry.node.id, entry]));
  const addIssue = (node, issue) => {
    const entry = byId.get(node?.id);
    if (entry) entry.issues.push(issue);
  };
  const globalIssues = [];

  activeEdges().forEach((edge) => {
    const from = nodeById(edge.from);
    const to = nodeById(edge.to);
    if (!from || !to) {
      globalIssues.push("A connection points to a missing block.");
      return;
    }
    if (edge.label && from.kind !== "data" && from.provides?.length && !from.provides.some((artifact) => artifactCompatible(edge.label, artifact))) {
      addIssue(from, `Connection label ${edge.label} is not documented as an output.`);
    }
    if (edge.label && to.kind !== "data" && to.requires?.length && !to.requires.some((artifact) => artifactCompatible(artifact, edge.label))) {
      addIssue(to, `Connection label ${edge.label} is not documented as an input.`);
    }
  });

  activeNodes().forEach((node) => {
    if (node.kind === "data") {
      if (!String(node.value || "").trim() && !activeEdges().some((edge) => edge.to === node.id)) {
        addIssue(node, "Initial data block has no value.");
      }
      return;
    }

    Object.entries(node.variables || {}).forEach(([name, value]) => {
      const cleanName = String(name || "").trim();
      if (!cleanName) {
        addIssue(node, "A selected variable has no name.");
        return;
      }
      if (!String(value || "").trim()) {
        addIssue(node, `Variable ${cleanName} has no value.`);
      }
    });

    Object.values(ensureInputValues(node)).forEach((input) => {
      const connected = incomingEdgesForArtifact(node, input.artifact).length > 0;
      const explicitValue = String(input.value || "").trim();
      if (explicitValue && !input.flag) {
        addIssue(node, `Input ${input.artifact} has a value but no parameter flag.`);
      }
      if (input.flag && !/^--?[A-Za-z0-9][A-Za-z0-9-]*$/.test(input.flag)) {
        addIssue(node, `Input ${input.artifact} has an invalid flag (${input.flag}).`);
      }
      if (!explicitValue && !connected && input.required) {
        addIssue(node, `Required input ${input.artifact} has no value or connection.`);
      }
    });

    (node.parameters || []).forEach((param) => {
      if (param.required && param.enabled && !String(param.value || param.defaultValue || "").trim()) {
        addIssue(node, `Required parameter ${param.name} has no value.`);
      }
      if (param.enabled && !String(param.name || "").trim()) {
        addIssue(node, "A selected parameter has no name.");
      }
      if (param.enabled && param.name && !/^--?[A-Za-z0-9][A-Za-z0-9-]*$/.test(param.name)) {
        addIssue(node, `Selected parameter ${param.name} is not a valid command-line flag.`);
      }
      if (param.enabled && !param.required && !parameterCanBeValueLess(param) && !String(param.value || "").trim()) {
        addIssue(node, `Selected parameter ${param.name} has no value.`);
      }
    });

    incompatibleSelectionIssues(node).forEach((issue) => addIssue(node, issue));
  });

  return [
    ...globalIssues.map((issue) => ({ node: null, title: "Workflow", kind: "workflow", issues: [issue] })),
    ...report
  ];
}

function parameterCanBeValueLess(param) {
  if (!param) return false;
  if (param.value != null && String(param.value).trim()) return false;
  const defaultValue = param.defaultValue == null ? "" : String(param.defaultValue).trim().toLowerCase();
  if (["true", "false"].includes(defaultValue)) return true;
  const docs = `${param.documentation || ""} ${param.description || ""}`.toLowerCase();
  return /\b(flag|toggle|skip|without|no-|only|report-only|overwrite)\b/.test(docs);
}

function selectedParameterSourceMap(node) {
  const sources = new Map();
  activeInputParameterNames(node).forEach((name) => {
    if (!sources.has(name)) sources.set(name, "connected input");
  });
  (node.parameters || []).forEach((param) => {
    if (!param.enabled) return;
    const label = param.name || "selected parameter";
    parameterNameCandidates(param).forEach((name) => {
      if (!sources.has(name)) sources.set(name, label);
    });
  });
  return sources;
}

function incompatibleSelectionIssues(node) {
  const sources = selectedParameterSourceMap(node);
  if (!sources.size) return [];
  const issues = [];
  (parameterRelationsForNode(node).incompatible || []).forEach((relation) => {
    const selected = relationParameters(relation).filter((name) => sources.has(name));
    if (selected.length < 2) return;
    const labels = [...new Set(selected.map((name) => sources.get(name) || name))];
    const evidence = relation?.evidence?.snippet ? ` (${compactText(relation.evidence.snippet, 110)})` : "";
    issues.push(`Incompatible parameters selected together: ${labels.join(", ")}.${evidence}`);
  });
  return [...new Set(issues)];
}

function showValidation() {
  const report = validateWorkflowReport();
  const issueCount = report.reduce((sum, entry) => sum + entry.issues.length, 0);
  const blockRows = report
    .filter((entry) => entry.node)
    .map((entry) => `
      <details class="validation-block" ${entry.issues.length ? "open" : ""}>
        <summary>
          <span>${escapeHtml(entry.title)}</span>
          <strong class="${entry.issues.length ? "validation-bad" : "validation-ok"}">${entry.issues.length ? `${entry.issues.length} issue${entry.issues.length === 1 ? "" : "s"}` : "OK"}</strong>
        </summary>
        ${entry.issues.length
          ? `<ul class="modal-list">${entry.issues.map((issue) => `<li>${escapeHtml(issue)}</li>`).join("")}</ul>`
          : `<p class="validation-ok">All selected variables and parameters are defined.</p>`}
      </details>
    `).join("");
  const globalRows = report
    .filter((entry) => !entry.node)
    .flatMap((entry) => entry.issues)
    .map((issue) => `<li>${escapeHtml(issue)}</li>`)
    .join("");

  if (!issueCount) {
    showModal("Workflow validation", `
      <p class="validation-ok">The current canvas has no broken links, missing explicitly required parameters, or empty selected variables.</p>
      <p>Official anvi'o pages expose many inputs as "can consume" relationships, so undocumented optional contexts are not forced.</p>
      <div class="validation-report">${blockRows}</div>
    `);
    return;
  }
  showModal("Workflow validation", `
    <p class="validation-bad">${issueCount} issue${issueCount === 1 ? "" : "s"} found.</p>
    ${globalRows ? `<ul class="modal-list">${globalRows}</ul>` : ""}
    <div class="validation-report">${blockRows}</div>
  `);
}

function workspacePayload() {
  return {
    schema: "org.anvio.workflow-builder.workspace",
    schemaVersion: "0.1.0",
    galaxyReady: true,
    createdWith: "anvio_workflow_site",
    documentationVersion: state.versionId,
    documentationSource: state.db.sources,
    settings: state.settings,
    tabs: state.tabs
  };
}

function saveWorkspace() {
  downloadBlob("anvio-workspace.json", JSON.stringify(workspacePayload(), null, 2), "application/json");
}

async function handleWorkspaceFile(event) {
  const file = event.target.files[0];
  event.target.value = "";
  if (!file) return;
  const text = await file.text();
  try {
    const payload = JSON.parse(text);
    if (!Array.isArray(payload.tabs)) throw new Error("Workspace JSON does not contain a tabs array.");
    withHistory(() => {
      state.tabs = payload.tabs;
      state.activeTabId = state.tabs[0]?.id || "";
      state.settings = {
        ...state.settings,
        ...(payload.settings || {})
      };
      if (payload.documentationVersion && state.db.versions.some((version) => version.id === payload.documentationVersion)) {
        state.versionId = payload.documentationVersion;
      }
      state.selectedNodeId = "";
      state.selectedEdgeId = "";
      applySettings();
    });
  } catch (error) {
    showModal("Load failed", `<p>${escapeHtml(error.message)}</p>`);
  }
}

async function handleWorkflowFile(event) {
  const file = event.target.files[0];
  event.target.value = "";
  if (!file) return;
  const text = await file.text();
  importWorkflowDraft = buildImportWorkflowDraft(file.name, text, importWorkflowDraft.selectedFormat || "auto");
  showImportWorkflowDialog();
}

function showImportWorkflowDialog() {
  showModal("Import workflow", `
    <div class="import-dialog">
      <div class="import-grid">
        <label class="field-stack">
          <span class="small-label">Workflow file</span>
          <input id="importWorkflowFileInput" type="file" accept=".sh,.bash,.smk,.snakefile,.txt,Snakefile">
        </label>
        <label class="field-stack">
          <span class="small-label">Format</span>
          <select id="importWorkflowFormat">
            <option value="auto">Auto-detect</option>
            <option value="bash">Bash</option>
            <option value="snakemake">Snakemake</option>
          </select>
        </label>
      </div>
      <div class="import-status" id="importWorkflowStatus">Choose a workflow file to preview it.</div>
      <pre class="import-preview" id="importWorkflowPreview">No file selected.</pre>
      <div class="import-command-list" id="importWorkflowCommands"></div>
      <div class="inspector-actions">
        <button class="secondary-button" type="button" data-action="confirm-import-workflow" id="confirmImportWorkflow" disabled>Import in new tab</button>
      </div>
    </div>
  `);
  const formatSelect = document.getElementById("importWorkflowFormat");
  if (formatSelect) formatSelect.value = importWorkflowDraft.selectedFormat || "auto";
  renderImportWorkflowDraft();
}

function openImportWorkflowDialog() {
  importWorkflowDraft = buildImportWorkflowDraft("", "", "auto");
  showImportWorkflowDialog();
}

async function handleModalChange(event) {
  const target = event.target;
  if (!(target instanceof Element)) return;

  if (target.id === "importWorkflowFileInput") {
    const file = target.files?.[0];
    if (!file) return;
    try {
      const text = await file.text();
      importWorkflowDraft = buildImportWorkflowDraft(file.name, text, importWorkflowDraft.selectedFormat || "auto");
      renderImportWorkflowDraft();
    } catch (error) {
      importWorkflowDraft = buildImportWorkflowDraft(file.name, "", importWorkflowDraft.selectedFormat || "auto");
      const status = document.getElementById("importWorkflowStatus");
      if (status) status.textContent = `Could not read file: ${error.message}`;
    }
  }

  if (target.id === "importWorkflowFormat") {
    importWorkflowDraft = buildImportWorkflowDraft(
      importWorkflowDraft.fileName,
      importWorkflowDraft.text,
      target.value || "auto"
    );
    renderImportWorkflowDraft();
  }
}

function buildImportWorkflowDraft(fileName, text, selectedFormat = "auto") {
  const detectedFormat = detectWorkflowFormat(fileName, text);
  const effectiveFormat = selectedFormat === "auto" ? detectedFormat : selectedFormat;
  const records = parseWorkflowCommandRecords(text, effectiveFormat);
  return {
    fileName,
    text,
    selectedFormat,
    detectedFormat,
    effectiveFormat,
    commands: records.map((record) => record.command),
    records,
    variables: detectWorkflowVariables(text, effectiveFormat)
  };
}

function detectWorkflowFormat(fileName, text) {
  const name = String(fileName || "").toLowerCase();
  if (name.endsWith(".smk") || name.endsWith(".snakefile") || name === "snakefile") return "snakemake";
  if (/\brule\s+[A-Za-z0-9_.-]+\s*:/.test(text) && /\bshell\s*:/.test(text)) return "snakemake";
  if (name.endsWith(".sh") || name.endsWith(".bash")) return "bash";
  if (/^#!.*\b(?:bash|sh)\b/m.test(text) || /^\s*(?:set\s+-[A-Za-z]*e|export\s+\w+=)/m.test(text)) return "bash";
  return "bash";
}

function formatLabel(format) {
  if (format === "snakemake") return "Snakemake";
  if (format === "bash") return "Bash";
  return "Unknown";
}

function detectWorkflowVariables(text, format = "auto") {
  const variables = {};
  const remember = (name, value = "") => {
    const cleanName = String(name || "").trim();
    if (!cleanName || /^\d+$/.test(cleanName)) return;
    if (!(cleanName in variables) || (!variables[cleanName] && value)) {
      variables[cleanName] = String(value || "").trim();
    }
  };

  text.split(/\r?\n/).forEach((line) => {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith("#")) return;

    const shellAssignment = trimmed.match(/^(?:export\s+)?([A-Za-z_][A-Za-z0-9_]*)=(.*)$/);
    if (shellAssignment) {
      remember(shellAssignment[1], cleanVariableValue(shellAssignment[2]));
    }

    const pythonAssignment = trimmed.match(/^([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.+)$/);
    if ((format === "snakemake" || format === "auto") && pythonAssignment && !/^(rule|shell|input|output|params|threads|resources)\b/.test(trimmed)) {
      remember(pythonAssignment[1], cleanVariableValue(pythonAssignment[2]));
    }

    const configFile = trimmed.match(/^configfile\s*:\s*(.+)$/);
    if (configFile) remember("configfile", cleanVariableValue(configFile[1]));
  });

  [...text.matchAll(/\$\{?([A-Za-z_][A-Za-z0-9_]*)[^A-Za-z0-9_}]?\}?/g)]
    .forEach((match) => remember(match[1]));

  [...text.matchAll(/\{(input|output|params|wildcards|resources|config)(?:\.([A-Za-z_][A-Za-z0-9_]*))?\}/g)]
    .forEach((match) => remember(match[2] ? `${match[1]}.${match[2]}` : match[1]));

  [...text.matchAll(/config\[['"]([^'"]+)['"]\]/g)]
    .forEach((match) => remember(`config.${match[1]}`));

  return Object.fromEntries(Object.entries(variables).sort(([a], [b]) => a.localeCompare(b)));
}

function cleanVariableValue(value) {
  return String(value || "")
    .replace(/\s+#.*$/, "")
    .replace(/^['"]|['"]$/g, "")
    .trim();
}

function renderImportWorkflowDraft() {
  const status = document.getElementById("importWorkflowStatus");
  const preview = document.getElementById("importWorkflowPreview");
  const commandsEl = document.getElementById("importWorkflowCommands");
  const confirmButton = document.getElementById("confirmImportWorkflow");
  if (!status || !preview || !commandsEl || !confirmButton) return;

  const draft = importWorkflowDraft;
  const variableEntries = Object.entries(draft.variables || {});
  if (!draft.text) {
    status.textContent = "Choose a workflow file to preview it.";
    preview.textContent = "No file selected.";
    commandsEl.innerHTML = "";
    confirmButton.disabled = true;
    return;
  }

  const detected = formatLabel(draft.detectedFormat);
  const effective = formatLabel(draft.effectiveFormat);
  const previewLimit = 8000;
  const clipped = draft.text.length > previewLimit;
  status.textContent = `${draft.fileName || "Workflow"} | detected: ${detected} | importing as: ${effective} | ${draft.commands.length} command${draft.commands.length === 1 ? "" : "s"} and ${variableEntries.length} variable${variableEntries.length === 1 ? "" : "s"} found.`;
  preview.textContent = `${draft.text.slice(0, previewLimit)}${clipped ? "\n\n... preview truncated ..." : ""}`;
  commandsEl.innerHTML = `
    ${draft.commands.length
      ? `<h3>Detected commands</h3><ul class="modal-list">${draft.commands.map((command) => `<li><code>${escapeHtml(command)}</code></li>`).join("")}</ul>`
      : `<p class="validation-bad">No known anvi'o or workflow external command was found in this file.</p>`}
    ${variableEntries.length
      ? `<h3>Detected variables</h3><ul class="modal-list">${variableEntries.map(([name, value]) => `<li><code>${escapeHtml(name)}</code>${value ? ` = <code>${escapeHtml(value)}</code>` : ""}</li>`).join("")}</ul>`
      : `<p>No shell or Snakemake variables were detected.</p>`}
  `;
  confirmButton.disabled = !draft.commands.length;
}

function confirmImportWorkflow() {
  const raw = importWorkflowDraft.records?.length
    ? importWorkflowDraft.records
    : (importWorkflowDraft.commands || []).map((command) => ({ command, raw: command, flags: [] }));
  if (!raw.length) {
    renderImportWorkflowDraft();
    return;
  }
  const seen = new Set();
  const records = raw.filter((record) => {
    if (!record.command || seen.has(record.command)) return false;
    seen.add(record.command);
    return true;
  });
  try {
    importWorkflowCommands(records, importWorkflowDraft.fileName || "Imported workflow", importWorkflowDraft.variables || {});
    if (typeof dom.modal.close === "function") dom.modal.close();
  } catch (error) {
    showModal("Import failed", `<p>${escapeHtml(error.message)}</p>`);
  }
}

function importWorkflowCommands(records, sourceName, variables = {}) {
  withHistory(() => {
    const programById = new Map(currentVersion().programs.map((program) => [program.id, program]));
    const importableRecords = records.filter((record) => record.command);
    if (!importableRecords.length) throw new Error("No command could be imported from this file.");

    const cleanName = String(sourceName || "Imported workflow")
      .replace(/\.[A-Za-z0-9_.-]+$/, "")
      .slice(0, 34) || "Imported workflow";
    const tab = createTab(cleanName);
    state.tabs.push(tab);
    state.activeTabId = tab.id;
    state.selectedNodeId = "";
    state.selectedEdgeId = "";
    state.connectStartId = "";

    const center = { x: WORLD_WIDTH / 2, y: WORLD_HEIGHT / 2 };
    const importedNodes = [];
    importableRecords.forEach((record, index) => {
      const item = programById.get(record.command) || externalWorkflowRecordToItem(record, sourceName);
      const node = spawnItem(item, center.x + index * 24, center.y + index * 24, false, false, false);
      if (node) {
        node.importOrder = index;
        applyImportedCommandRecord(node, record, variables);
        importedNodes.push(node);
      }
    });
    connectImportedWorkflowNodes(importedNodes);
    addImportedInitialInputDataBlocks(importedNodes, center);
    applyOptimizedPlacement(center);
  });
  requestAnimationFrame(centerCanvasViewport);
}

function externalWorkflowRecordToItem(record, sourceName = "Imported workflow") {
  const parameters = (record.flags || []).map(({ flag, value }) => ({
    name: flag,
    value: value || "",
    defaultValue: null,
    required: false,
    enabled: true,
    source: "imported-external-workflow",
    documentation: "Detected while importing an external workflow command."
  }));
  return {
    id: `external:${record.command}`,
    kind: "program",
    title: record.command,
    command: record.command,
    description: "External command detected in an imported workflow.",
    requires: [],
    provides: [],
    parameters,
    examples: [record.raw || record.command],
    docsUrl: GITHUB_URL,
    documentation: `External command imported from ${sourceName}. It is not part of the local anvi'o runtime database, so its parameters are fully editable.`,
    externalFromAnvio: true,
    externalSource: sourceName
  };
}

function addImportedInitialInputDataBlocks(importedNodes, center) {
  const providedByArtifact = new Map();
  importedNodes.forEach((node) => {
    (node.provides || []).forEach((artifact) => {
      if (!providedByArtifact.has(artifact)) providedByArtifact.set(artifact, node);
    });
  });
  const initialInputs = new Map();

  importedNodes.forEach((node) => {
    Object.values(ensureInputValues(node)).forEach((input) => {
      const value = String(input.value || "").trim();
      if (!value) return;
      const producer = providedByArtifact.get(input.artifact);
      if (producer && producer.id !== node.id) {
        ensureEdge(producer.id, node.id, input.artifact);
        return;
      }
      const key = `${input.artifact}\u0000${value}`;
      if (!initialInputs.has(key)) {
        initialInputs.set(key, { artifact: input.artifact, value, consumers: [] });
      }
      initialInputs.get(key).consumers.push(node);
    });
  });

  [...initialInputs.values()].forEach((input, index) => {
    const artifactId = input.artifact;
    const item = findDataBlock(artifactId);
    if (!item) return;
    const dataNode = createNodeFromItem(item, center.x - 360, center.y + index * 140);
    dataNode.importOrder = Math.max(0, Math.min(...input.consumers.map((node) => node.importOrder ?? index)) - 0.5);
    dataNode.value = input.value;
    normalizeNodeSize(dataNode);
    currentTab().nodes.push(dataNode);
    input.consumers.forEach((consumer) => ensureEdge(dataNode.id, consumer.id, artifactId));
  });
}

function connectImportedWorkflowNodes(importedNodes) {
  importedNodes.forEach((consumer, consumerIndex) => {
    (consumer.requires || []).forEach((artifactId) => {
      const producer = importedNodes
        .slice(0, consumerIndex)
        .filter((candidate) => nodeProvidesArtifact(candidate, artifactId))
        .at(-1);
      if (producer) ensureEdge(producer.id, consumer.id, artifactId);
    });
  });
}

function applyImportedCommandRecord(node, record, workflowVariables = {}) {
  const recordVariables = detectWorkflowVariables(record.raw || "", importWorkflowDraft.effectiveFormat || "auto");
  node.variables = {
    ...(node.variables || {}),
    ...clone(workflowVariables),
    ...recordVariables
  };

  const inputs = ensureInputValues(node);
  (record.flags || []).forEach(({ flag, value }) => {
    let matchingInput = importedFlagInputMatch(node, inputs, flag);
    if (!matchingInput) {
      const artifactId = importedFlagArtifactId(flag);
      if (artifactId) {
        if (!node.requires.includes(artifactId)) node.requires.push(artifactId);
        ensureInputValues(node);
        matchingInput = node.inputValues[artifactId];
      }
    }
    if (matchingInput) {
      matchingInput.value = value;
      return;
    }

    let param = (node.parameters || []).find((candidate) => candidate.name === flag);
    if (!param) {
      param = {
        name: flag,
        value: "",
        defaultValue: null,
        required: false,
        enabled: true,
        source: "imported-workflow",
        documentation: "Detected while importing a workflow file."
      };
      node.parameters.push(param);
    }
    param.value = value;
    param.enabled = true;
  });

  enforceParameterCompatibility(node);
  normalizeNodeSize(node);
}

function importedFlagArtifactId(flag) {
  const normalizedFlag = normalizedDataKey(flag);
  if (!normalizedFlag) return "";
  const dataBlocks = currentVersion().dataBlocks || [];
  const exact = dataBlocks.find((item) => normalizedDataKey(item.id) === normalizedFlag);
  if (exact) return exact.id;
  const aliases = {
    c: "contigs-db",
    "contigs-db": "contigs-db",
    i: "input",
    input: "input",
    f: "fasta",
    fasta: "fasta",
    p: "profile-db",
    "profile-db": "profile-db",
    pan: "pan-db",
    "pan-db": "pan-db",
    g: "genomes-storage-db",
    "genomes-storage": "genomes-storage-db",
    "genomes-storage-db": "genomes-storage-db",
    "external-genomes": "external-genomes",
    "internal-genomes": "internal-genomes",
    "external-gene-calls": "external-gene-calls"
  };
  const aliased = aliases[normalizedFlag];
  if (!aliased || aliased === "input" || aliased === "fasta") return "";
  return dataBlocks.some((item) => item.id === aliased) ? aliased : "";
}

function importedFlagInputMatch(node, inputs, flag) {
  const inputList = Object.values(inputs || {});
  const exact = inputList.find((input) => input.flag === flag);
  if (exact) return exact;

  const aliases = {
    "-c": ["contigs-db"],
    "--contigs-db": ["contigs-db"],
    "-i": ["bam-file", "raw-bam-file", "single-profile-db", "profile-db"],
    "--input": ["bam-file", "raw-bam-file", "single-profile-db", "profile-db"],
    "-f": ["contigs-fasta", "fasta", "fasta-txt"],
    "--fasta": ["contigs-fasta", "fasta", "fasta-txt"],
    "-p": ["profile-db", "single-profile-db", "pan-db"],
    "--profile-db": ["profile-db", "single-profile-db"],
    "-P": ["pan-db"],
    "--pan-db": ["pan-db"],
    "-g": ["genomes-storage-db"],
    "--genomes-storage": ["genomes-storage-db"],
    "--genomes-storage-db": ["genomes-storage-db"],
    "--external-genomes": ["external-genomes"],
    "--internal-genomes": ["internal-genomes"]
  };
  const artifactHints = aliases[flag] || [];
  const aliasMatch = artifactHints
    .map((artifact) => inputList.find((input) => input.artifact === artifact))
    .find(Boolean);
  if (aliasMatch) return aliasMatch;

  const flagWords = normalizedDataKey(flag).split("-").filter((word) => word.length > 1);
  return inputList.find((input) => {
    const artifactWords = normalizedDataKey(input.artifact).split("-");
    return flagWords.some((word) => artifactWords.includes(word));
  }) || null;
}

function parseWorkflowCommandRecords(text, format = "auto") {
  const documented = new Set(currentVersion().programs.map((program) => program.id));
  const knownCommands = workflowImportCommandSet();
  const normalized = String(text || "").replace(/\\\r?\n/g, " ");
  const records = [];
  const matcher = /\b([A-Za-z0-9_.-]+)\b([^\n;]*)/g;
  let match;

  while ((match = matcher.exec(normalized))) {
    const command = match[1];
    if (!documented.has(command) && !knownCommands.has(command)) continue;
    const argsText = match[2] || "";
    records.push({
      command,
      raw: `${command}${argsText}`,
      flags: parseCommandFlags(argsText)
    });
  }

  return records;
}

function workflowImportCommandSet() {
  const commands = new Set(currentVersion().programs.map((program) => program.id));
  if (Array.isArray(window.ANVIO_PRELOADED_WORKFLOWS)) {
    window.ANVIO_PRELOADED_WORKFLOWS.forEach((workflow) => {
      (workflow.steps || []).forEach((step) => {
        if (step.command) commands.add(step.command);
      });
    });
  }
  return commands;
}

function parseCommandFlags(argsText) {
  const tokens = tokenizeCommandArgs(argsText);
  const flags = [];
  for (let index = 0; index < tokens.length; index += 1) {
    const token = tokens[index];
    if (!/^--?[A-Za-z0-9][A-Za-z0-9-]*$/.test(token)) continue;
    const next = tokens[index + 1] || "";
    const hasValue = next && !next.startsWith("-");
    flags.push({ flag: token, value: hasValue ? next : "" });
    if (hasValue) index += 1;
  }
  return flags;
}

function tokenizeCommandArgs(argsText) {
  return [...String(argsText || "").matchAll(/"([^"]*)"|'([^']*)'|`([^`]*)`|(\S+)/g)]
    .map((match) => match[1] ?? match[2] ?? match[3] ?? match[4] ?? "")
    .map((token) => token.replace(/^['"`]|['"`]$/g, ""))
    .map((token) => token.replace(/[),]+$/, ""))
    .filter(Boolean);
}

function parseWorkflowCommands(text, format = "auto") {
  return parseWorkflowCommandRecords(text, format)
    .map((record) => record.command)
    .filter((command, index, array) => array.indexOf(command) === index);
}

function openExportDialog() {
  const tab = currentTab();
  const defaultTitle = tab?.name || "anvio workflow";
  showModal("Export current canvas", `
    <div class="export-dialog">
      <label class="field-stack" for="exportTitleInput">
        <span>Title</span>
        <input id="exportTitleInput" type="text" value="${escapeHtml(defaultTitle)}" placeholder="Workflow title">
      </label>
      <div class="export-format-grid" aria-label="Export formats">
        <button class="secondary-button export-format-button" type="button" data-action="export-format" data-format="bash">
          <svg viewBox="0 0 24 24" class="menu-icon" aria-hidden="true"><path d="M4 5h16v14H4z"></path><path d="m7 9 3 3-3 3"></path><path d="M12 15h5"></path></svg>
          <span>Bash script</span>
        </button>
        <button class="secondary-button export-format-button" type="button" data-action="export-format" data-format="snakemake">
          <svg viewBox="0 0 24 24" class="menu-icon" aria-hidden="true"><path d="M5 4h14v16H5z"></path><path d="M8 8h8"></path><path d="M8 12h8"></path><path d="M8 16h5"></path></svg>
          <span>Snakemake</span>
        </button>
        <button class="secondary-button export-format-button" type="button" data-action="export-format" data-format="svg">
          <svg viewBox="0 0 24 24" class="menu-icon" aria-hidden="true"><path d="M4 5h16v14H4z"></path><path d="M8 16 6 12l2-4"></path><path d="m16 8 2 4-2 4"></path><path d="m14 7-4 10"></path></svg>
          <span>Graph SVG</span>
        </button>
        <button class="secondary-button export-format-button" type="button" data-action="export-format" data-format="png">
          <svg viewBox="0 0 24 24" class="menu-icon" aria-hidden="true"><rect x="4" y="5" width="16" height="14" rx="2"></rect><circle cx="9" cy="10" r="1.5"></circle><path d="m7 17 4-4 3 3 2-2 3 3"></path></svg>
          <span>Graph PNG</span>
        </button>
      </div>
    </div>
  `);
}

function exportSelectedFormat(event) {
  const button = event.target.closest("[data-format]");
  const format = button?.dataset.format || "";
  const title = document.getElementById("exportTitleInput")?.value.trim() || currentTab()?.name || "anvio workflow";
  const exporters = {
    bash: () => exportBash(title),
    snakemake: () => exportSnakemake(title),
    svg: () => exportSvg(title),
    png: () => exportPng(title)
  };
  exporters[format]?.();
  if (typeof dom.modal.close === "function") dom.modal.close();
}

function exportFilename(title, extension) {
  const base = safeName(title || currentTab()?.name || "anvio-workflow").replace(/_/g, "-") || "anvio-workflow";
  return `${base}.${extension}`;
}

function exportBash(title = currentTab()?.name || "anvio workflow") {
  const lines = [
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "",
    `# ${title}`,
    `# Generated from ${currentVersion().label} local anvi'o runtime metadata.`,
    ""
  ];
  orderedExecutableNodes().forEach((node) => {
    lines.push(`# ${node.title}`);
    if (node.description) lines.push(`# ${compactText(node.description, 110)}`);
    lines.push(buildCommand(node));
    lines.push("");
  });
  downloadBlob(exportFilename(title, "sh"), lines.join("\n"), "text/x-shellscript");
}

function exportSnakemake(title = currentTab()?.name || "anvio workflow") {
  const rules = orderedExecutableNodes().map((node, index) => {
    const incoming = activeEdges().filter((edge) => edge.to === node.id);
    const outgoing = activeEdges().filter((edge) => edge.from === node.id);
    const inputLines = incoming.length
      ? incoming.map((edge) => edgeToSnakemakeLine(edge, "from")).join(",\n")
      : "        ";
    const outputLines = outgoing.length
      ? outgoing.map((edge) => edgeToSnakemakeLine(edge, "to")).join(",\n")
      : "        ";
    return [
      `rule ${safeName(node.title)}_${index + 1}:`,
      "    input:",
      inputLines,
      "    output:",
      outputLines,
      "    shell:",
      `        ${JSON.stringify(buildCommand(node))}`
    ].join("\n");
  });
  downloadBlob(exportFilename(title, "smk"), `# ${title}\n# Generated from ${currentVersion().label} local anvi'o runtime metadata.\n\n${rules.join("\n\n")}`, "text/plain");
}

function edgeToSnakemakeLine(edge, side) {
  const node = nodeById(edge[side]);
  const artifact = edge.label || node?.refId || "artifact";
  const value = node?.kind === "data" ? (node.value || node.refId) : artifact;
  return `        ${safeName(artifact)}=${JSON.stringify(value)}`;
}

function orderedExecutableNodes() {
  const executable = activeNodes().filter((node) => node.kind !== "data");
  const dependencies = new Map(executable.map((node) => [node.id, new Set()]));
  activeEdges().forEach((edge) => {
    const from = nodeById(edge.from);
    const to = nodeById(edge.to);
    if (!from || !to) return;
    if (from.kind !== "data" && to.kind !== "data") dependencies.get(to.id)?.add(from.id);
    if (from.kind !== "data" && to.kind === "data") {
      activeEdges()
        .filter((next) => next.from === to.id)
        .map((next) => nodeById(next.to))
        .filter((nextNode) => nextNode && nextNode.kind !== "data")
        .forEach((nextNode) => dependencies.get(nextNode.id)?.add(from.id));
    }
  });

  const result = [];
  const remaining = new Map([...dependencies].map(([id, deps]) => [id, new Set(deps)]));
  while (remaining.size) {
    const ready = [...remaining].find(([, deps]) => deps.size === 0);
    if (!ready) break;
    const [id] = ready;
    result.push(executable.find((node) => node.id === id));
    remaining.delete(id);
    remaining.forEach((deps) => deps.delete(id));
  }
  const unresolved = executable.filter((node) => !result.includes(node));
  return [...result, ...unresolved].filter(Boolean);
}

function exportSvg(title = currentTab()?.name || "anvio workflow") {
  downloadBlob(exportFilename(title, "svg"), graphSvgMarkup({ title, cropped: true }), "image/svg+xml");
}

function exportPng(title = currentTab()?.name || "anvio workflow") {
  const bounds = graphBounds();
  const svg = graphSvgMarkup({ title, cropped: true, bounds });
  const blob = new Blob([svg], { type: "image/svg+xml" });
  const url = URL.createObjectURL(blob);
  const image = new Image();
  image.onload = () => {
    const canvas = document.createElement("canvas");
    canvas.width = Math.ceil(bounds.width);
    canvas.height = Math.ceil(bounds.height);
    const context = canvas.getContext("2d");
    context.fillStyle = getComputedStyle(document.body).getPropertyValue("--bg").trim() || "#ffffff";
    context.fillRect(0, 0, canvas.width, canvas.height);
    context.drawImage(image, 0, 0);
    URL.revokeObjectURL(url);
    canvas.toBlob((pngBlob) => {
      downloadBlob(exportFilename(title, "png"), pngBlob, "image/png");
    });
  };
  image.src = url;
}

function graphBounds(padding = 48) {
  const nodes = activeNodes();
  const edges = activeEdges();
  if (!nodes.length) return { x: 0, y: 0, width: 900, height: 600 };
  const points = [];
  nodes.forEach((node) => {
    points.push({ x: node.x, y: node.y });
    points.push({ x: node.x + node.w, y: node.y + node.h });
  });
  edges.forEach((edge) => routeEdge(edge).forEach((point) => points.push(point)));
  const minX = Math.max(0, Math.min(...points.map((point) => point.x)) - padding);
  const minY = Math.max(0, Math.min(...points.map((point) => point.y)) - padding);
  const maxX = Math.min(WORLD_WIDTH, Math.max(...points.map((point) => point.x)) + padding);
  const maxY = Math.min(WORLD_HEIGHT, Math.max(...points.map((point) => point.y)) + padding);
  return {
    x: minX,
    y: minY,
    width: Math.max(240, maxX - minX),
    height: Math.max(180, maxY - minY)
  };
}

function graphSvgMarkup(options = {}) {
  const bounds = options.bounds || (options.cropped ? graphBounds() : { x: 0, y: 0, width: WORLD_WIDTH, height: WORLD_HEIGHT });
  const nodes = activeNodes();
  const edges = activeEdges();
  const styles = getComputedStyle(document.body);
  const bg = styles.getPropertyValue("--bg").trim() || "#f5f1e8";
  const panel = styles.getPropertyValue("--panel").trim() || "#fffdf8";
  const ink = styles.getPropertyValue("--ink").trim() || "#20252b";
  const accent = styles.getPropertyValue("--accent-strong").trim() || "#145f68";
  const line = styles.getPropertyValue("--line-strong").trim() || "#aeb9c5";
  const blocks = nodes.map((node) => `
    <g>
      <rect x="${node.x}" y="${node.y}" width="${node.w}" height="${node.h}" rx="8" fill="${panel}" stroke="${line}" stroke-width="2"></rect>
      <text x="${node.x + 10}" y="${node.y + 24}" fill="${ink}" font-family="Arial" font-size="14" font-weight="700">${escapeHtml(compactText(node.title, 28))}</text>
      <text x="${node.x + 10}" y="${node.y + 48}" fill="${ink}" font-family="Arial" font-size="11">${escapeHtml(compactText(node.kind === "data" ? node.artifactType : node.command, 34))}</text>
    </g>
  `).join("");
  const arrows = edges.map((edge) => {
    const d = pointsToPath(routeEdge(edge));
    return `<path d="${d}" fill="none" stroke="${accent}" stroke-width="2.2" marker-end="url(#arrowhead-export)"></path>`;
  }).join("");
  const title = options.title ? `<title>${escapeHtml(options.title)}</title>` : "";
  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${Math.ceil(bounds.width)}" height="${Math.ceil(bounds.height)}" viewBox="${bounds.x} ${bounds.y} ${bounds.width} ${bounds.height}">
  ${title}
  <defs>
    <marker id="arrowhead-export" markerWidth="10" markerHeight="8" refX="9" refY="4" orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L10,4 L0,8 z" fill="${accent}"></path>
    </marker>
  </defs>
  <rect x="${bounds.x}" y="${bounds.y}" width="${bounds.width}" height="${bounds.height}" fill="${bg}"></rect>
  ${arrows}
  ${blocks}
</svg>`;
}

function downloadBlob(filename, content, type) {
  const blob = content instanceof Blob ? content : new Blob([content], { type });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = filename;
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function showHelp() {
  showModal("Guide", `
    <p>This builder uses the official anvi'o help pages as a searchable block database. Search descriptions on the left, then click a program, workflow, or data block to place it on the active canvas.</p>
    <ul class="modal-list">
      <li><strong>File:</strong> save or load the entire workspace as JSON; import bash or Snakemake files to detect documented anvi'o commands.</li>
      <li><strong>Edition:</strong> undo, redo, adjust font size, and toggle night mode.</li>
      <li><strong>Canvas:</strong> drag blocks, resize from the lower-right corner, drag from one block port to another block to create arrows, drag empty canvas space to pan, and use zoom or optimize placement from the toolbar.</li>
      <li><strong>Inspector:</strong> edit titles, descriptions, variables, documented parameters, custom parameters, and open official documentation.</li>
      <li><strong>Export:</strong> produce bash, Snakemake, SVG, or PNG from the active canvas.</li>
    </ul>
    <p>Example: search for <code>contigs database</code>, add <code>anvi-gen-contigs-database</code>, set documented or custom parameters in the inspector, then validate and export.</p>
  `);
}

function showCitation() {
  showModal("Citation", `
    <p>Cite anvi'o using the official citation guidance from the anvi'o project. The app links directly to the canonical citation page below.</p>
    <p><a href="https://anvio.org/cite/" target="_blank" rel="noopener">https://anvio.org/cite/</a></p>
    <p>This local builder database is generated from the installed anvi'o runtime and <code>ARCHITECTURE.md</code>.</p>
  `);
}

function openOfficialDocPopup(event) {
  const button = event.target.closest("[data-doc-url]");
  const url = button?.dataset.docUrl;
  if (!url) return;
  const node = activeNodes().find((candidate) => candidate.id === state.selectedNodeId);
  const libraryItem = node ? findLibraryItem(node.kind, node.refId) : null;
  const title = libraryItem?.title || node?.title || "Official anvi'o documentation";
  const description = libraryItem?.description || node?.description || "Official anvi'o documentation page.";
  const documentation = libraryItem?.documentation || node?.documentation || description;

  showModal("Official anvi'o documentation", `
    <div class="official-doc-popup integrated-doc">
      <div class="official-doc-source">
        <span>Official source</span>
        <a href="${escapeHtml(url)}" target="_blank" rel="noopener">${escapeHtml(url)}</a>
      </div>
      <article class="official-doc-page">
        <h3>${escapeHtml(title)}</h3>
        <p>${escapeHtml(description)}</p>
        <pre class="official-doc-content">${escapeHtml(documentation)}</pre>
      </article>
    </div>
  `);
}

function showModal(title, html) {
  dom.modalContent.innerHTML = `<h2>${escapeHtml(title)}</h2>${html}`;
  if (typeof dom.modal.showModal === "function") dom.modal.showModal();
  else alert(`${title}\n\n${dom.modalContent.textContent}`);
}

init();
