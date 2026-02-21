use std::collections::BTreeMap;

/// Records a summary of transcript interactions for debugging.
///
/// Port of C++ `TranscriptManifest` from transcript_manifest.hpp.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TranscriptManifest {
    manifest: BTreeMap<usize, RoundData>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RoundData {
    challenge_labels: Vec<String>,
    entries: Vec<(String, usize)>,
}

impl Default for TranscriptManifest {
    fn default() -> Self {
        Self::new()
    }
}

impl TranscriptManifest {
    pub fn new() -> Self {
        Self {
            manifest: BTreeMap::new(),
        }
    }

    /// Add a single challenge label to the manifest for the given round.
    pub fn add_challenge(&mut self, round: usize, label: &str) {
        self.manifest
            .entry(round)
            .or_insert_with(|| RoundData {
                challenge_labels: Vec::new(),
                entries: Vec::new(),
            })
            .challenge_labels
            .push(label.to_string());
    }

    /// Add an entry (element label + field count) to the manifest for the given round.
    pub fn add_entry(&mut self, round: usize, label: &str, element_size: usize) {
        self.manifest
            .entry(round)
            .or_insert_with(|| RoundData {
                challenge_labels: Vec::new(),
                entries: Vec::new(),
            })
            .entries
            .push((label.to_string(), element_size));
    }

    pub fn size(&self) -> usize {
        self.manifest.len()
    }

    pub fn print(&self) {
        for (round, data) in &self.manifest {
            println!("Round: {round}");
            for label in &data.challenge_labels {
                println!("\tchallenge: {label}");
            }
            for (label, size) in &data.entries {
                println!("\telement ({size}): {label}");
            }
        }
    }
}
