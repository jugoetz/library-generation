import os
from typing import List, Dict

from labware.plates import Plate, Plate96, Plate384, Plate384Echo


class TransferStep:
    def __init__(
        self,
        source_id: str,
        source_well: str,
        destination_id: str,
        destination_well: str,
        volume: int,
    ):
        """
        :param source_id: source plate identifier
        :param source_well: source well for the transfer, e.g. 'A1'
        :param destination_id: destination plate identifier
        :param destination_well: destination well for the transfer, e.g. 'A1'
        :param volume: transfer volume in nL
        """
        self.source_id = source_id
        self.source_well = source_well
        self.destination_id = destination_id
        self.destination_well = destination_well
        self.volume = int(volume)  # in nL

    def __repr__(self):
        return f"TransferStep[{self.source_id}_{self.source_well} --> {self.destination_id}_{self.destination_well}, {self.volume} nL]"


class Transfer:
    def __init__(self, transfer_steps: List[TransferStep]):
        self.transfer_steps = transfer_steps
        self.step_count = len(transfer_steps)

    def __repr__(self):
        s = "Transfer:\n"
        for i in self.transfer_steps:
            s += f"{i}\n"
        s += f"Total of {self.step_count} steps"
        return s

    def simulate(
        self,
        source_plates: Dict[str, Plate],
        destination_well_number: int,
        save_plate=False,
        save_dir="",
    ):
        """
        Simulate the outcome of the Transfer, given a source plate.

        Args:
            source_plates (dict): Dictionary of the source plates of the form {"name": labware.plates.Plate}
            destination_well_number (int): Well count of the destination plate. Can be 96 or 384.
            save_plate (bool): Whether to save the resulting destination plate. Default False.
            save_dir (str): Path to save the destination plate.
        """

        destination_plates = {}
        for transfer_step in self.transfer_steps:
            # check if source plate referenced in Echo file is in source_plates
            if transfer_step.source_id not in source_plates.keys():
                raise KeyError("Source plate not found")
            else:
                source_plate = source_plates[transfer_step.source_id]

            # check if well referenced in Echo file exists
            if transfer_step.source_well not in source_plate.wells():
                raise KeyError("Invalid source well")
            else:
                source_cmps, source_vol = source_plate.well(transfer_step.source_well)

            # remove from source well if volume allows
            try:
                source_plate.consume_well(
                    transfer_step.source_well, transfer_step.volume
                )
            except ValueError:
                raise ValueError(
                    f"Cannot retrieve {transfer_step.volume} nL from well {transfer_step.source_id}_{transfer_step.source_well} filled with {source_vol} nL"
                )

            # add destination plate to dict if it was not encountered before
            if transfer_step.destination_id not in destination_plates.keys():
                if destination_well_number == 96:
                    destination_plates[transfer_step.destination_id] = Plate96(
                        max_vol=9999999999, dead_vol=0
                    )
                elif destination_well_number == 384:
                    destination_plates[transfer_step.destination_id] = Plate384(
                        max_vol=9999999999, dead_vol=0
                    )
                else:
                    raise NotImplementedError("Can only use 96 or 384 well plates")

            # add to destination plate
            try:
                destination_plates[transfer_step.destination_id].fill_well(
                    transfer_step.destination_well, source_cmps, transfer_step.volume
                )
            except Exception:
                ValueError(
                    f"Invalid transfer {transfer_step} for given destination plate"
                )

            # save destination plates if requested
            if save_plate is True:
                if save_dir == "":
                    raise ValueError("save_plate set to True but no save_dir given.")
                else:
                    # save plates to file
                    for name, plate in destination_plates.items():
                        plate.to_csv(
                            os.path.join(save_dir, f"plate_layout_{name}.csv"),
                            save_volumes=True,
                        )

            else:
                if save_dir != "":
                    raise ValueError(
                        "save_dir given but save_plate set to False. Set save_plate to True if you want to save the destination plate."
                    )

        return destination_plates

    def check_correct(
        self, source_plates: Dict[str, Plate], destination_plates: Dict[str, Plate]
    ):
        pass


if __name__ == "__main__":
    """
    This routine checks whether given Echo transfer files generate the correct target plates from a given source plate.
    """

    import csv

    def import_echo_transfer_file(file):
        with open(file, "r") as file:
            reader = csv.reader(file)
            transfer_steps = []
            next(reader, None)  # skip the header
            for line in reader:
                transfer_steps.append(
                    TransferStep(line[0], line[1], line[2], line[3], int(line[4]))
                )
        return transfer_steps

    transfer_steps = []
    transfer_steps += import_echo_transfer_file("../data/plates/exp25/step1.csv")
    transfer_steps += import_echo_transfer_file("../data/plates/exp25/step2.csv")
    transfer = Transfer(transfer_steps)
    print(transfer)
    source_plate = Plate384(
        max_vol=65000, dead_vol=13000
    )  # note we set the dead volume to 13 uL not 15 as in plate specs
    source_plate.from_csv(
        "/Users/julian/PycharmProjects/library-generation/data/plates/exp25/source_plate_layout.csv"
    )
    print(source_plate)
    planned_destination_plates = {
        "Synthesis1": Plate384Echo(),
        "Synthesis2": Plate384Echo(),
        "Synthesis3": Plate384Echo(),
        "Synthesis4": Plate384Echo(),
        "Synthesis5": Plate384Echo(),
        "Synthesis6": Plate384Echo(),
    }

    for i in range(6):
        planned_destination_plates[f"Synthesis{i + 1}"].from_csv(
            f"/Users/julian/PycharmProjects/library-generation/data/plates/exp25/plate_layout_plate{i + 1}.csv"
        )
        print((planned_destination_plates[f"Synthesis{i + 1}"]))

    destination_plates = transfer.simulate(
        {"Source1": source_plate},
        384,
        save_plate=True,
        save_dir="/Users/julian/Desktop",
    )

    for key, plate in destination_plates.items():
        planned_plate = planned_destination_plates[key]
        print(plate)
        for well, well_planned in zip(
            plate.iterate_wells(), planned_plate.iterate_wells()
        ):
            well_planned[1].insert(2, "X")
            if well != well_planned:
                print(f"{well}, {well_planned}")
