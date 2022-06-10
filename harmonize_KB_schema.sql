CREATE TABLE IF NOT EXISTS "harmonization_version_control" (
    pgs_id VARCHAR(25) NOT NULL,
    version INT NOT NULL,
    genebuild INT NOT NULL,
    date TEXT NOT NULL,
    source TEXT NOT NULL,
    variant_number INT NOT NULL,
    matched_number INT NOT NULL,
    unmapped_number INT NOT NULL,
    comment TEXT NOT NULL,
    PRIMARY KEY (pgs_id, version, genebuild)
);